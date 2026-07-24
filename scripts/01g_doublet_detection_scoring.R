#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
})

option_list <- list(
  make_option(
    "--input",
    type = "character",
    default = NULL,
    help = "Input QC-filtered Seurat RDS object"
  ),
  make_option(
    "--config",
    type = "character",
    default = NULL,
    help = "YAML doublet detection configuration file"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_doublet_scored.rds",
    help = "Output Seurat RDS object annotated with doublet scores"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

message_step <- function(label) {
  message("\n========== ", label, " ==========")
}

check_file <- function(path, label) {
  if (is.null(path) || !nzchar(path)) {
    stop("Missing path for ", label)
  }

  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
}

get_config_value <- function(config, path, default = NULL) {
  parts <- strsplit(path, "\\.")[[1]]
  value <- config

  for (part in parts) {
    if (is.null(value[[part]])) {
      return(default)
    }

    value <- value[[part]]
  }

  value
}

get_doubletfinder_function <- function(base_name) {
  candidates <- c(
    paste0(base_name, "_v3"),
    base_name
  )

  for (fn in candidates) {
    if (
      exists(
        fn,
        where = asNamespace("DoubletFinder"),
        inherits = FALSE
      )
    ) {
      return(get(fn, envir = asNamespace("DoubletFinder")))
    }

    if (exists(fn, mode = "function")) {
      return(get(fn, mode = "function"))
    }
  }

  stop("DoubletFinder function not found: ", base_name)
}

safe_load_doubletfinder <- function() {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'DoubletFinder' is not installed in the container. ",
        "Install it in the Docker image before running this module."
      )
    )
  }

  suppressPackageStartupMessages({
    library(DoubletFinder)
  })

  paramSweep_fun <<- get_doubletfinder_function("paramSweep")
  summarizeSweep_fun <<- get_doubletfinder_function("summarizeSweep")
  findPK_fun <<- get_doubletfinder_function("find.pK")
  modelHomotypic_fun <<- get_doubletfinder_function("modelHomotypic")
  doubletFinder_fun <<- get_doubletfinder_function("doubletFinder")
}

call_with_supported_args <- function(fun, args) {
  fn_formals <- names(formals(fun))
  args <- args[names(args) %in% fn_formals]
  do.call(fun, args)
}

get_counts_layers <- function(seu, assay = DefaultAssay(seu)) {
  assay_obj <- seu[[assay]]
  available_layers <- Layers(assay_obj)
  counts_layers <- grep("^counts", available_layers, value = TRUE)

  if (length(counts_layers) == 0) {
    stop(
      paste0(
        "No counts-like layers found in assay ",
        assay,
        ". Available layers: ",
        paste(available_layers, collapse = ", ")
      )
    )
  }

  counts_layers
}

get_sample_counts_matrix <- function(
    seu,
    sample_id,
    sample_cells,
    assay = DefaultAssay(seu)
) {
  assay_obj <- seu[[assay]]
  counts_layers <- get_counts_layers(seu, assay = assay)

  direct_layers <- grep(
    pattern = paste0("^counts\\.", sample_id, "(\\.|$)"),
    x = counts_layers,
    value = TRUE
  )

  candidate_layers <- if (length(direct_layers) > 0) {
    direct_layers
  } else {
    counts_layers
  }

  best_layer <- NULL
  best_overlap <- 0

  for (layer_name in candidate_layers) {
    mat_tmp <- LayerData(
      object = assay_obj,
      layer = layer_name
    )

    overlap <- length(
      intersect(sample_cells, colnames(mat_tmp))
    )

    if (overlap > best_overlap) {
      best_overlap <- overlap
      best_layer <- layer_name
    }

    rm(mat_tmp)
    gc()
  }

  if (is.null(best_layer) || best_overlap == 0) {
    stop("No matching counts layer found for sample: ", sample_id)
  }

  mat <- LayerData(
    object = assay_obj,
    layer = best_layer
  )

  mat <- as(mat, "dgCMatrix")
  common_cells <- intersect(sample_cells, colnames(mat))

  if (length(common_cells) == 0) {
    stop(
      "No shared cells between metadata and counts layer for sample: ",
      sample_id
    )
  }

  mat[, common_cells, drop = FALSE]
}

safe_scale_and_pca <- function(obj, npcs, nfeatures) {
  obj <- NormalizeData(
    object = obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  variable_features <- VariableFeatures(obj)

  obj <- ScaleData(
    object = obj,
    features = variable_features,
    verbose = FALSE
  )

  RunPCA(
    object = obj,
    features = variable_features,
    npcs = npcs,
    verbose = FALSE
  )
}

select_best_pk <- function(sweep_stats) {
  bcmvn_raw <- findPK_fun(sweep_stats)

  if (is.null(bcmvn_raw)) {
    stop(
      "DoubletFinder find.pK returned NULL instead of a BCmvn table."
    )
  }

  bcmvn <- as.data.frame(
    bcmvn_raw,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(bcmvn) == 0) {
    stop("DoubletFinder returned an empty BCmvn table.")
  }

  colnames(bcmvn) <- make.unique(colnames(bcmvn))

  pk_index <- grep(
    "^pK($|\\.)",
    colnames(bcmvn)
  )

  bc_index <- grep(
    "^BCmetric($|\\.)",
    colnames(bcmvn)
  )

  if (length(pk_index) == 0) {
    stop(
      "Column pK not found in BCmvn table. Available columns: ",
      paste(colnames(bcmvn), collapse = ", ")
    )
  }

  if (length(bc_index) == 0) {
    stop(
      "Column BCmetric not found in BCmvn table. Available columns: ",
      paste(colnames(bcmvn), collapse = ", ")
    )
  }

  pK_numeric <- suppressWarnings(
    as.numeric(
      as.character(
        bcmvn[[pk_index[[1]]]]
      )
    )
  )

  BCmetric_numeric <- suppressWarnings(
    as.numeric(
      as.character(
        bcmvn[[bc_index[[1]]]]
      )
    )
  )

  valid_indices <- which(
    !is.na(pK_numeric) &
      is.finite(pK_numeric) &
      !is.na(BCmetric_numeric) &
      is.finite(BCmetric_numeric)
  )

  if (length(valid_indices) == 0) {
    stop(
      "No valid numeric pK and BCmetric combinations found."
    )
  }

  best_index <- valid_indices[
    which.max(
      BCmetric_numeric[valid_indices]
    )
  ]

  best_pk <- as.numeric(
    pK_numeric[[best_index]]
  )

  if (
    length(best_pk) != 1 ||
      is.na(best_pk) ||
      !is.finite(best_pk)
  ) {
    stop("Selected pK is not a valid numeric scalar.")
  }

  bcmvn$pK_numeric <- pK_numeric
  bcmvn$BCmetric_numeric <- BCmetric_numeric

  message(
    "Selected pK: ",
    best_pk,
    " | BCmetric: ",
    BCmetric_numeric[[best_index]]
  )

  list(
    best_pk = best_pk,
    bcmvn = bcmvn
  )
}

get_doubletfinder_columns <- function(meta_before, meta_after) {
  new_cols <- setdiff(colnames(meta_after), colnames(meta_before))

  classification_col <- grep(
    pattern = "^DF.classifications",
    x = new_cols,
    value = TRUE
  )

  pann_col <- grep(
    pattern = "^pANN",
    x = new_cols,
    value = TRUE
  )

  if (length(classification_col) == 0) {
    stop("DoubletFinder classification column not found.")
  }

  if (length(pann_col) == 0) {
    stop("DoubletFinder pANN column not found.")
  }

  list(
    classification_col = classification_col[length(classification_col)],
    pann_col = pann_col[length(pann_col)]
  )
}

skipped_sample_result <- function(
    sample_id,
    cell_ids,
    expected_doublet_rate,
    pN,
    status
) {
  metadata <- data.frame(
    cell_id = cell_ids,
    sample = sample_id,
    doublet_score = NA_real_,
    doublet_class = "NotRun",
    is_doublet = NA,
    doublet_pK = NA_real_,
    doublet_nExp_raw = NA_integer_,
    doublet_nExp_adjusted = NA_integer_,
    doublet_detection_status = "not_run",
    doublet_decision_status = "unreviewed",
    doublet_reviewed_class = NA_character_,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    sample = sample_id,
    cells_scored = length(cell_ids),
    expected_doublet_rate = expected_doublet_rate,
    pN = pN,
    pK = NA_real_,
    nExp_raw = NA_integer_,
    homotypic_prop = NA_real_,
    nExp_adjusted = NA_integer_,
    proposed_singlets = NA_integer_,
    proposed_doublets = NA_integer_,
    doublet_detection_status = status,
    doublet_decision_status = "unreviewed",
    approval_status = "proposed",
    stringsAsFactors = FALSE
  )

  list(
    metadata = metadata,
    summary = summary,
    bcmvn = NULL
  )
}

run_doubletfinder_one_sample <- function(
    sample_id,
    seu,
    meta,
    pcs,
    expected_doublet_rate,
    pN,
    min_cells,
    nfeatures,
    reuse_pANN = NULL
) {
  message_step(paste0("Scoring sample: ", sample_id))

  sample_cells <- rownames(meta)[meta$sample == sample_id]

  if (length(sample_cells) < min_cells) {
    message(
      "Skipping sample ",
      sample_id,
      ": cells below min_cells threshold."
    )

    return(
      skipped_sample_result(
        sample_id = sample_id,
        cell_ids = sample_cells,
        expected_doublet_rate = expected_doublet_rate,
        pN = pN,
        status = "skipped_low_cell_number"
      )
    )
  }

  sample_counts <- get_sample_counts_matrix(
    seu = seu,
    sample_id = sample_id,
    sample_cells = sample_cells,
    assay = DefaultAssay(seu)
  )

  sample_meta <- meta[colnames(sample_counts), , drop = FALSE]

  obj <- CreateSeuratObject(
    counts = sample_counts,
    project = sample_id,
    meta.data = sample_meta
  )

  rm(sample_counts)
  gc()

  DefaultAssay(obj) <- "RNA"

  max_pcs <- min(pcs, ncol(obj) - 1, nrow(obj) - 1)

  if (max_pcs < 5) {
    cell_ids <- colnames(obj)
    rm(obj)
    gc()

    return(
      skipped_sample_result(
        sample_id = sample_id,
        cell_ids = cell_ids,
        expected_doublet_rate = expected_doublet_rate,
        pN = pN,
        status = "skipped_too_few_pcs"
      )
    )
  }

  obj <- safe_scale_and_pca(
    obj = obj,
    npcs = max_pcs,
    nfeatures = nfeatures
  )

  pcs_use <- seq_len(max_pcs)

  sweep_res <- call_with_supported_args(
    fun = paramSweep_fun,
    args = list(
      seu = obj,
      PCs = pcs_use,
      sct = FALSE,
      pN = pN
    )
  )

  sweep_stats <- summarizeSweep_fun(
    sweep_res,
    GT = FALSE
  )

  pk_res <- select_best_pk(sweep_stats)
  best_pk <- pk_res$best_pk

  if (
    !is.numeric(best_pk) ||
      length(best_pk) != 1 ||
      is.na(best_pk) ||
      !is.finite(best_pk)
  ) {
    stop(
      "Invalid pK selected for sample ",
      sample_id,
      ": expected one finite numeric value."
    )
  }

  best_pk <- as.numeric(best_pk)

  bcmvn <- pk_res$bcmvn
  bcmvn$sample <- sample_id

  nExp_raw <- max(
    1,
    round(ncol(obj) * expected_doublet_rate)
  )

  homotypic_prop <- if (
    "seurat_clusters" %in% colnames(obj@meta.data)
  ) {
    tryCatch(
      modelHomotypic_fun(obj$seurat_clusters),
      error = function(e) 0
    )
  } else {
    0
  }

  nExp_adjusted <- max(
    1,
    round(nExp_raw * (1 - homotypic_prop))
  )

  meta_before <- obj@meta.data

  obj <- call_with_supported_args(
    fun = doubletFinder_fun,
    args = list(
      seu = obj,
      PCs = pcs_use,
      pN = pN,
      pK = best_pk,
      nExp = nExp_adjusted,
      reuse.pANN = reuse_pANN,
      sct = FALSE
    )
  )

  df_cols <- get_doubletfinder_columns(
    meta_before = meta_before,
    meta_after = obj@meta.data
  )

  doublet_class_raw <- as.character(
    obj@meta.data[[df_cols$classification_col]]
  )

  doublet_score <- as.numeric(
    obj@meta.data[[df_cols$pann_col]]
  )

  doublet_class <- ifelse(
    doublet_class_raw == "Doublet",
    "Doublet",
    "Singlet"
  )

  is_doublet <- doublet_class == "Doublet"

  sample_metadata <- data.frame(
    cell_id = colnames(obj),
    sample = sample_id,
    doublet_score = doublet_score,
    doublet_class = doublet_class,
    is_doublet = is_doublet,
    doublet_pK = best_pk,
    doublet_nExp_raw = nExp_raw,
    doublet_nExp_adjusted = nExp_adjusted,
    doublet_detection_status = "completed",
    doublet_decision_status = "unreviewed",
    doublet_reviewed_class = NA_character_,
    stringsAsFactors = FALSE
  )

  sample_summary <- data.frame(
    sample = sample_id,
    cells_scored = ncol(obj),
    expected_doublet_rate = expected_doublet_rate,
    pN = pN,
    pK = best_pk,
    nExp_raw = nExp_raw,
    homotypic_prop = homotypic_prop,
    nExp_adjusted = nExp_adjusted,
    proposed_singlets = sum(!is_doublet, na.rm = TRUE),
    proposed_doublets = sum(is_doublet, na.rm = TRUE),
    doublet_detection_status = "completed",
    doublet_decision_status = "unreviewed",
    approval_status = "proposed",
    stringsAsFactors = FALSE
  )

  rm(obj)
  gc()

  list(
    metadata = sample_metadata,
    summary = sample_summary,
    bcmvn = bcmvn
  )
}

check_file(opt$input, "input Seurat object")
check_file(opt$config, "doublet configuration")

message_step("Loading DoubletFinder")
safe_load_doubletfinder()

message_step("Loading configuration")
config <- yaml::read_yaml(opt$config)

expected_doublet_rate <- get_config_value(
  config,
  "doublet.expected_rate",
  0.075
)

pN <- get_config_value(config, "doublet.pN", 0.25)
pcs <- get_config_value(config, "doublet.pcs", 15)
min_cells <- get_config_value(config, "doublet.min_cells", 500)
nfeatures <- get_config_value(config, "doublet.nfeatures", 1500)

message_step("Loading QC-filtered Seurat object")
seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

if (ncol(seu) == 0) {
  stop("Input Seurat object contains no cells")
}

DefaultAssay(seu) <- "RNA"

if (!"sample" %in% colnames(seu@meta.data)) {
  if ("ShortName" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$ShortName
  } else if ("sample_id" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$sample_id
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$orig.ident
  } else {
    stop(
      paste(
        "No sample metadata column found.",
        "Expected one of: sample, ShortName, sample_id, orig.ident"
      )
    )
  }
}

if (!"Group" %in% colnames(seu@meta.data)) {
  seu$Group <- NA_character_
}

meta <- seu@meta.data
meta$cell_id <- rownames(meta)
samples <- sort(unique(as.character(meta$sample)))

message("Samples detected: ", paste(samples, collapse = ", "))

message_step("Running per-sample doublet scoring")
set.seed(1234)

results <- lapply(
  samples,
  run_doubletfinder_one_sample,
  seu = seu,
  meta = meta,
  pcs = pcs,
  expected_doublet_rate = expected_doublet_rate,
  pN = pN,
  min_cells = min_cells,
  nfeatures = nfeatures,
  reuse_pANN = NULL
)

doublet_metadata <- dplyr::bind_rows(
  lapply(results, `[[`, "metadata")
)

doublet_summary_by_sample <- dplyr::bind_rows(
  lapply(results, `[[`, "summary")
)

bcmvn_list <- lapply(results, `[[`, "bcmvn")
bcmvn_list <- bcmvn_list[!sapply(bcmvn_list, is.null)]

if (length(bcmvn_list) > 0) {
  bcmvn_table <- dplyr::bind_rows(bcmvn_list)
} else {
  bcmvn_table <- data.frame(
    sample = character(),
    pK = numeric(),
    BCmetric = numeric()
  )
}

message_step("Adding doublet annotations without filtering")

rownames(doublet_metadata) <- doublet_metadata$cell_id
doublet_metadata <- doublet_metadata[colnames(seu), , drop = FALSE]

if (!identical(rownames(doublet_metadata), colnames(seu))) {
  stop("Doublet metadata does not match Seurat cell identity/order")
}

seu$doublet_score <- doublet_metadata$doublet_score
seu$doublet_class <- doublet_metadata$doublet_class
seu$is_doublet <- doublet_metadata$is_doublet
seu$doublet_detection_status <-
  doublet_metadata$doublet_detection_status
seu$doublet_decision_status <-
  doublet_metadata$doublet_decision_status
seu$doublet_reviewed_class <-
  doublet_metadata$doublet_reviewed_class
seu$doublet_pK <- doublet_metadata$doublet_pK
seu$doublet_nExp_raw <- doublet_metadata$doublet_nExp_raw
seu$doublet_nExp_adjusted <- doublet_metadata$doublet_nExp_adjusted

cells_input <- ncol(seu)
features_input <- nrow(seu)
cells_output <- ncol(seu)
features_output <- nrow(seu)
cells_removed <- 0
proposed_doublets <- sum(
  seu$is_doublet,
  na.rm = TRUE
)

overall_summary <- data.frame(
  metric = c(
    "cells_input",
    "cells_output",
    "features_input",
    "features_output",
    "cells_removed",
    "samples",
    "proposed_doublets",
    "proposed_doublet_fraction",
    "expected_doublet_rate",
    "pN",
    "pcs",
    "min_cells",
    "nfeatures",
    "approval_status"
  ),
  value = c(
    cells_input,
    cells_output,
    features_input,
    features_output,
    cells_removed,
    length(samples),
    proposed_doublets,
    proposed_doublets / cells_input,
    expected_doublet_rate,
    pN,
    pcs,
    min_cells,
    nfeatures,
    "proposed"
  )
)

message_step("Saving scoring outputs")

saveRDS(seu, opt$output)

readr::write_csv(
  overall_summary,
  "01g_doublet_detection_summary.csv"
)

readr::write_csv(
  doublet_summary_by_sample,
  "01g_doublet_detection_by_sample.csv"
)

readr::write_csv(
  doublet_metadata,
  "01g_doublet_detection_cell_scores.csv"
)

readr::write_csv(
  bcmvn_table,
  "01g_doublet_detection_pk_sweep.csv"
)

sink("01g_sessionInfo.txt")
print(sessionInfo())
sink()

message("Doublet-scored object: ", opt$output)
message("No cells were removed by this module.")
message_step("01g_doublet_detection_scoring finished successfully")
