#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
  library(ggplot2)
})

# -----------------------------
# Arguments
# -----------------------------

option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "Input QC-filtered Seurat RDS object"
  ),
  make_option(
    c("--config"),
    type = "character",
    default = NULL,
    help = "YAML configuration file"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = "seu_singlets.rds",
    help = "Output Seurat RDS object after doublet removal"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Missing --input")
}

if (is.null(opt$config)) {
  stop("Missing --config")
}

# -----------------------------
# Helper functions
# -----------------------------

message_step <- function(x) {
  message("\n========== ", x, " ==========")
}

check_file <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(paste0("File does not exist - ", label, ": ", path))
  }
}

get_config_value <- function(config, path, default = NULL) {
  parts <- strsplit(path, "\\.")[[1]]
  value <- config

  for (p in parts) {
    if (is.null(value[[p]])) {
      return(default)
    }
    value <- value[[p]]
  }

  value
}

get_doubletfinder_function <- function(base_name) {
  candidates <- c(
    paste0(base_name, "_v3"),
    base_name
  )

  for (fn in candidates) {
    if (exists(fn, where = asNamespace("DoubletFinder"), inherits = FALSE)) {
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

  message(
    "Available layers in assay ",
    assay,
    ": ",
    paste(available_layers, collapse = ", ")
  )

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

  message("Counts layers detected: ", paste(counts_layers, collapse = ", "))

  counts_layers
}

get_sample_counts_matrix <- function(seu, sample_id, sample_cells, assay = DefaultAssay(seu)) {

  assay_obj <- seu[[assay]]
  counts_layers <- get_counts_layers(seu, assay = assay)

  direct_layers <- grep(
    pattern = paste0("^counts\\.", sample_id, "(\\.|$)"),
    x = counts_layers,
    value = TRUE
  )

  candidate_layers <- if (length(direct_layers) > 0) direct_layers else counts_layers

  best_layer <- NULL
  best_overlap <- 0

  for (layer_name in candidate_layers) {

    message("Testing layer for sample ", sample_id, ": ", layer_name)

    mat_tmp <- LayerData(
      object = assay_obj,
      layer = layer_name
    )

    overlap <- length(intersect(sample_cells, colnames(mat_tmp)))

    message("Cell overlap with ", sample_id, ": ", overlap)

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

  message("Using layer for sample ", sample_id, ": ", best_layer)

  mat <- LayerData(
    object = assay_obj,
    layer = best_layer
  )

  mat <- as(mat, "dgCMatrix")

  common_cells <- intersect(sample_cells, colnames(mat))

  if (length(common_cells) == 0) {
    stop("No shared cells between metadata and counts layer for sample: ", sample_id)
  }

  mat <- mat[, common_cells, drop = FALSE]

  message(
    "Sample ",
    sample_id,
    " matrix: ",
    nrow(mat),
    " genes x ",
    ncol(mat),
    " cells"
  )

  mat
}

safe_scale_and_pca <- function(obj, npcs = 20, nfeatures = 1500) {

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

  obj <- RunPCA(
    object = obj,
    features = variable_features,
    npcs = npcs,
    verbose = FALSE
  )

  obj
}

select_best_pk <- function(sweep_stats) {

  bcmvn <- findPK_fun(sweep_stats)

  bcmvn$pK_numeric <- as.numeric(as.character(bcmvn$pK))
  bcmvn$BCmetric <- as.numeric(as.character(bcmvn$BCmetric))

  bcmvn <- bcmvn %>%
    dplyr::filter(!is.na(pK_numeric), !is.na(BCmetric))

  if (nrow(bcmvn) == 0) {
    stop("No valid pK values found by DoubletFinder.")
  }

  best_row <- bcmvn[which.max(bcmvn$BCmetric), ]

  list(
    best_pk = best_row$pK_numeric,
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

run_doubletfinder_one_sample <- function(
  sample_id,
  seu,
  meta,
  pcs,
  expected_doublet_rate,
  pN,
  min_cells,
  nfeatures,
  reuse_pANN = FALSE
) {

  message_step(paste0("Processing sample: ", sample_id))

  sample_cells <- rownames(meta)[meta$sample == sample_id]

  message("Cells in sample before doublet detection: ", length(sample_cells))

  if (length(sample_cells) < min_cells) {

    message(
      "Sample ",
      sample_id,
      " has fewer than ",
      min_cells,
      " cells. Doublet detection skipped."
    )

    skipped_metadata <- data.frame(
      cell_id = sample_cells,
      sample = sample_id,
      doublet_score = NA_real_,
      doublet_class = "NotRun",
      is_doublet = FALSE,
      pK = NA_real_,
      nExp_raw = NA_integer_,
      nExp_adjusted = NA_integer_,
      stringsAsFactors = FALSE
    )

    skipped_summary <- data.frame(
      sample = sample_id,
      cells_before = length(sample_cells),
      expected_doublet_rate = expected_doublet_rate,
      pN = pN,
      pK = NA_real_,
      nExp_raw = NA_integer_,
      homotypic_prop = NA_real_,
      nExp_adjusted = NA_integer_,
      singlets = length(sample_cells),
      doublets = 0,
      status = "skipped_low_cell_number",
      stringsAsFactors = FALSE
    )

    return(
      list(
        metadata = skipped_metadata,
        summary = skipped_summary,
        bcmvn = NULL
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

    message("Too few cells/features for PCA. Doublet detection skipped for sample: ", sample_id)

    skipped_metadata <- data.frame(
      cell_id = colnames(obj),
      sample = sample_id,
      doublet_score = NA_real_,
      doublet_class = "NotRun",
      is_doublet = FALSE,
      pK = NA_real_,
      nExp_raw = NA_integer_,
      nExp_adjusted = NA_integer_,
      stringsAsFactors = FALSE
    )

    skipped_summary <- data.frame(
      sample = sample_id,
      cells_before = ncol(obj),
      expected_doublet_rate = expected_doublet_rate,
      pN = pN,
      pK = NA_real_,
      nExp_raw = NA_integer_,
      homotypic_prop = NA_real_,
      nExp_adjusted = NA_integer_,
      singlets = ncol(obj),
      doublets = 0,
      status = "skipped_too_few_pcs",
      stringsAsFactors = FALSE
    )

    rm(obj)
    gc()

    return(
      list(
        metadata = skipped_metadata,
        summary = skipped_summary,
        bcmvn = NULL
      )
    )
  }

  obj <- safe_scale_and_pca(
    obj = obj,
    npcs = max_pcs,
    nfeatures = nfeatures
  )

  pcs_use <- seq_len(max_pcs)

  message("Running DoubletFinder parameter sweep")

  param_sweep_args <- list(
    seu = obj,
    PCs = pcs_use,
    sct = FALSE,
    pN = pN
  )

  sweep_res <- call_with_supported_args(
    fun = paramSweep_fun,
    args = param_sweep_args
  )

  sweep_stats <- summarizeSweep_fun(
    sweep_res,
    GT = FALSE
  )

  pk_res <- select_best_pk(sweep_stats)

  best_pk <- pk_res$best_pk
  bcmvn <- pk_res$bcmvn
  bcmvn$sample <- sample_id

  message("Selected pK: ", best_pk)

  nExp_raw <- round(ncol(obj) * expected_doublet_rate)

  if (nExp_raw < 1) {
    nExp_raw <- 1
  }

  # Homotypic correction only if clusters are already present.
  # At this stage, clusters usually are not available, so use 0 to avoid under-calling doublets.
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    homotypic_prop <- tryCatch(
      {
        modelHomotypic_fun(obj$seurat_clusters)
      },
      error = function(e) {
        message("modelHomotypic failed. Using homotypic proportion = 0.")
        0
      }
    )
  } else {
    homotypic_prop <- 0
  }

  nExp_adjusted <- round(nExp_raw * (1 - homotypic_prop))

  if (nExp_adjusted < 1) {
    nExp_adjusted <- 1
  }

  message("nExp raw: ", nExp_raw)
  message("homotypic proportion: ", homotypic_prop)
  message("nExp adjusted: ", nExp_adjusted)

  meta_before <- obj@meta.data

  doubletfinder_args <- list(
    seu = obj,
    PCs = pcs_use,
    pN = pN,
    pK = best_pk,
    nExp = nExp_adjusted,
    reuse.pANN = reuse_pANN,
    sct = FALSE
  )

  obj <- call_with_supported_args(
    fun = doubletFinder_fun,
    args = doubletfinder_args
  )

  df_cols <- get_doubletfinder_columns(
    meta_before = meta_before,
    meta_after = obj@meta.data
  )

  doublet_class_raw <- as.character(obj@meta.data[[df_cols$classification_col]])
  doublet_score <- as.numeric(obj@meta.data[[df_cols$pann_col]])

  doublet_class <- ifelse(doublet_class_raw == "Doublet", "Doublet", "Singlet")
  is_doublet <- doublet_class == "Doublet"

  sample_metadata <- data.frame(
    cell_id = colnames(obj),
    sample = sample_id,
    doublet_score = doublet_score,
    doublet_class = doublet_class,
    is_doublet = is_doublet,
    pK = best_pk,
    nExp_raw = nExp_raw,
    nExp_adjusted = nExp_adjusted,
    stringsAsFactors = FALSE
  )

  singlets <- sum(!is_doublet, na.rm = TRUE)
  doublets <- sum(is_doublet, na.rm = TRUE)

  sample_summary <- data.frame(
    sample = sample_id,
    cells_before = ncol(obj),
    expected_doublet_rate = expected_doublet_rate,
    pN = pN,
    pK = best_pk,
    nExp_raw = nExp_raw,
    homotypic_prop = homotypic_prop,
    nExp_adjusted = nExp_adjusted,
    singlets = singlets,
    doublets = doublets,
    status = "completed",
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

save_doublet_plots <- function(metadata, summary_table) {

  p1 <- ggplot(metadata, aes(x = sample, fill = doublet_class)) +
    geom_bar() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Doublet classification by sample",
      x = "Sample",
      y = "Cells",
      fill = "Class"
    )

  ggsave(
    filename = "01c_doublet_class_by_sample.png",
    plot = p1,
    width = 8,
    height = 5,
    dpi = 150
  )

  if (sum(!is.na(metadata$doublet_score)) > 0) {
    p2 <- ggplot(metadata, aes(x = sample, y = doublet_score)) +
      geom_boxplot(outlier.size = 0.3) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Doublet score by sample",
        x = "Sample",
        y = "Doublet score"
      )

    ggsave(
      filename = "01c_doublet_score_by_sample.png",
      plot = p2,
      width = 8,
      height = 5,
      dpi = 150
    )
  }

  p3 <- ggplot(summary_table, aes(x = sample, y = doublets)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Detected doublets by sample",
      x = "Sample",
      y = "Detected doublets"
    )

  ggsave(
    filename = "01c_doublet_count_by_sample.png",
    plot = p3,
    width = 8,
    height = 5,
    dpi = 150
  )
}

# -----------------------------
# Load packages
# -----------------------------

message_step("Loading DoubletFinder")

safe_load_doubletfinder()

# -----------------------------
# Load configuration
# -----------------------------

message_step("Loading configuration")

check_file(opt$config, "configuration file")
config <- yaml::read_yaml(opt$config)

expected_doublet_rate <- get_config_value(config, "doublet.expected_rate", 0.075)
pN <- get_config_value(config, "doublet.pN", 0.25)
pcs <- get_config_value(config, "doublet.pcs", 15)
min_cells <- get_config_value(config, "doublet.min_cells", 500)
nfeatures <- get_config_value(config, "doublet.nfeatures", 1500)

message("expected_doublet_rate: ", expected_doublet_rate)
message("pN: ", pN)
message("pcs: ", pcs)
message("min_cells: ", min_cells)
message("nfeatures: ", nfeatures)

# -----------------------------
# Load object
# -----------------------------

message_step("Loading QC-filtered Seurat object")

check_file(opt$input, "input Seurat object")

seu <- readRDS(opt$input)

message("Input cells: ", ncol(seu))
message("Input genes: ", nrow(seu))
message("Default assay: ", DefaultAssay(seu))

DefaultAssay(seu) <- "RNA"

# -----------------------------
# Ensure sample metadata
# -----------------------------

message_step("Checking metadata")

if (!"sample" %in% colnames(seu@meta.data)) {
  if ("ShortName" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$ShortName
  } else if ("sample_id" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$sample_id
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$orig.ident
  } else {
    stop("No sample metadata column found. Expected one of: sample, ShortName, sample_id, orig.ident")
  }
}

if (!"Group" %in% colnames(seu@meta.data)) {
  seu$Group <- NA_character_
}

meta <- seu@meta.data
meta$cell_id <- rownames(meta)

# -----------------------------
# Prepare metadata only
# -----------------------------

message_step("Preparing per-sample doublet detection")

message("Cells available for doublet detection: ", nrow(meta))
message("Genes available in Seurat object: ", nrow(seu))

samples <- sort(unique(as.character(meta$sample)))

message("Samples detected: ", paste(samples, collapse = ", "))

# -----------------------------
# Run DoubletFinder per sample
# -----------------------------

message_step("Running DoubletFinder per sample")

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
  reuse_pANN = FALSE
)

doublet_metadata <- dplyr::bind_rows(lapply(results, `[[`, "metadata"))
doublet_summary_by_sample <- dplyr::bind_rows(lapply(results, `[[`, "summary"))

bcmvn_list <- lapply(results, `[[`, "bcmvn")
bcmvn_list <- bcmvn_list[!sapply(bcmvn_list, is.null)]

if (length(bcmvn_list) > 0) {
  bcmvn_table <- dplyr::bind_rows(bcmvn_list)
  readr::write_csv(bcmvn_table, "01c_doublet_pk_sweep.csv")
}

# -----------------------------
# Add metadata to original object
# -----------------------------

message_step("Adding doublet metadata to Seurat object")

rownames(doublet_metadata) <- doublet_metadata$cell_id

doublet_metadata <- doublet_metadata[colnames(seu), , drop = FALSE]

seu$doublet_score <- doublet_metadata$doublet_score
seu$doublet_class <- doublet_metadata$doublet_class
seu$is_doublet <- doublet_metadata$is_doublet
seu$doublet_pK <- doublet_metadata$pK
seu$doublet_nExp_raw <- doublet_metadata$nExp_raw
seu$doublet_nExp_adjusted <- doublet_metadata$nExp_adjusted

# -----------------------------
# Remove doublets
# -----------------------------

message_step("Removing doublets")

cells_before <- ncol(seu)

singlet_cells <- colnames(seu)[
  is.na(seu$is_doublet) |
    seu$is_doublet == FALSE
]

seu_singlets <- subset(
  x = seu,
  cells = singlet_cells
)

cells_after <- ncol(seu_singlets)
cells_removed <- cells_before - cells_after

message("Cells before doublet removal: ", cells_before)
message("Cells after doublet removal: ", cells_after)
message("Doublets removed: ", cells_removed)

# -----------------------------
# Summary tables
# -----------------------------

message_step("Saving summaries")

overall_summary <- data.frame(
  metric = c(
    "cells_before_doublet_removal",
    "cells_after_doublet_removal",
    "doublets_removed",
    "expected_doublet_rate",
    "pN",
    "pcs",
    "min_cells",
    "nfeatures"
  ),
  value = c(
    cells_before,
    cells_after,
    cells_removed,
    expected_doublet_rate,
    pN,
    pcs,
    min_cells,
    nfeatures
  )
)

readr::write_csv(overall_summary, "01c_doublet_summary.csv")
readr::write_csv(doublet_summary_by_sample, "01c_doublet_by_sample.csv")
readr::write_csv(doublet_metadata, "01c_doublet_metadata.csv")

save_doublet_plots(
  metadata = doublet_metadata,
  summary_table = doublet_summary_by_sample
)

# -----------------------------
# Save object
# -----------------------------

message_step("Saving singlet object")

saveRDS(seu_singlets, opt$output)

sink("01c_sessionInfo.txt")
print(sessionInfo())
sink()

message("Singlet object saved in: ", opt$output)
message("Doublet summary saved in: 01c_doublet_summary.csv")
message("Doublet by sample saved in: 01c_doublet_by_sample.csv")
message("Doublet metadata saved in: 01c_doublet_metadata.csv")
message("Session info saved in: 01c_sessionInfo.txt")

message_step("01c_doublet_detection finished successfully")