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
  ),
  make_option(
    "--cores",
    type = "integer",
    default = NULL,
    help = paste(
      "CPU allocation supplied by the workflow.",
      "Overrides doublet.sweep_cores when provided."
    )
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

message_step <- function(label) {
  message("\n========== ", label, " ==========")
}

time_step <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start

  list(
    value = value,
    seconds = as.numeric(elapsed)
  )
}

message_timed <- function(sample_id, label, seconds) {
  message(
    "[",
    sample_id,
    "] ",
    label,
    " finished in ",
    round(seconds, 2),
    " s"
  )
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
  t_norm <- time_step(
    NormalizeData(
      object = obj,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  )
  obj <- t_norm$value

  t_hvg <- time_step(
    FindVariableFeatures(
      object = obj,
      selection.method = "vst",
      nfeatures = nfeatures,
      verbose = FALSE
    )
  )
  obj <- t_hvg$value

  variable_features <- VariableFeatures(obj)

  t_scale <- time_step(
    ScaleData(
      object = obj,
      features = variable_features,
      verbose = FALSE
    )
  )
  obj <- t_scale$value

  t_pca <- time_step(
    RunPCA(
      object = obj,
      features = variable_features,
      npcs = npcs,
      verbose = FALSE
    )
  )

  list(
    object = t_pca$value,
    variable_features = variable_features,
    timings = list(
      normalize_seconds = t_norm$seconds,
      find_variable_features_seconds = t_hvg$seconds,
      scale_data_seconds = t_scale$seconds,
      run_pca_seconds = t_pca$seconds
    )
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
    best_bcmetric = as.numeric(BCmetric_numeric[[best_index]]),
    best_index = best_index,
    bcmvn = bcmvn
  )
}

resolve_sweep_size <- function(total_cells, min_cells, fraction, max_cells) {
  target <- floor(total_cells * fraction)
  target <- max(min_cells, target)
  target <- min(max_cells, target, total_cells)
  as.integer(target)
}

allocate_stratified_quota <- function(stratum_sizes, target_cells) {
  stratum_names <- names(stratum_sizes)
  stratum_sizes <- as.integer(stratum_sizes)
  names(stratum_sizes) <- stratum_names
  if (is.null(names(stratum_sizes)) || any(!nzchar(names(stratum_sizes)))) {
    stop("Stratum sizes must preserve non-empty names.")
  }

  exact_quota <- target_cells * stratum_sizes / sum(stratum_sizes)
  names(exact_quota) <- names(stratum_sizes)
  quota <- floor(exact_quota)
  names(quota) <- names(stratum_sizes)
  quota <- pmax(1L, quota)
  names(quota) <- names(stratum_sizes)
  quota <- pmin(quota, stratum_sizes)
  names(quota) <- names(stratum_sizes)

  while (sum(quota) > target_cells) {
    candidates <- names(quota)[quota > 1L]
    if (length(candidates) == 0) break
    chosen <- candidates[[which.max(quota[candidates] - exact_quota[candidates])]]
    quota[[chosen]] <- quota[[chosen]] - 1L
  }

  while (sum(quota) < target_cells) {
    capacity <- stratum_sizes - quota
    names(capacity) <- names(quota)
    candidates <- names(capacity)[capacity > 0L]
    if (length(candidates) == 0) break
    residual <- exact_quota[candidates] - quota[candidates]
    chosen <- candidates[[which.max(residual)]]
    quota[[chosen]] <- quota[[chosen]] + 1L
  }

  if (sum(quota) != target_cells) {
    stop("Unable to allocate the requested stratified sample size.")
  }

  quota
}

pca_stratified_cells <- function(
    obj,
    target_cells,
    pcs_use,
    seed,
    max_strata = 20L
) {
  all_cells <- colnames(obj)

  if (target_cells >= length(all_cells)) {
    return(all_cells)
  }

  embeddings <- Embeddings(obj, reduction = "pca")
  dimensions <- intersect(pcs_use, seq_len(ncol(embeddings)))

  if (length(dimensions) == 0) {
    set.seed(seed)
    return(sample(all_cells, target_cells, replace = FALSE))
  }

  coordinates <- embeddings[all_cells, dimensions, drop = FALSE]
  unique_coordinates <- nrow(unique(as.data.frame(coordinates)))
  if (unique_coordinates < 2L) {
    set.seed(seed)
    return(sample(all_cells, target_cells, replace = FALSE))
  }

  strata_count <- min(
    as.integer(max_strata),
    max(2L, as.integer(floor(sqrt(target_cells / 25)))),
    nrow(coordinates),
    unique_coordinates
  )

  set.seed(seed)
  clusters <- stats::kmeans(
    coordinates,
    centers = strata_count,
    iter.max = 50,
    nstart = 1
  )$cluster

  strata <- as.character(clusters)
  names(strata) <- rownames(coordinates)
  stratum_sizes <- table(strata)
  quota <- allocate_stratified_quota(stratum_sizes, target_cells)

  set.seed(seed)
  selected <- unlist(
    lapply(names(quota), function(stratum) {
      candidates <- names(strata)[strata == stratum]
      sample(candidates, quota[[stratum]], replace = FALSE)
    }),
    use.names = FALSE
  )

  if (length(selected) != target_cells) {
    stop(
      "PCA-stratified sampling returned ",
      length(selected),
      " cells instead of ",
      target_cells,
      "."
    )
  }

  selected
}

run_pk_sweep <- function(
    obj,
    pcs_use,
    pN,
    sweep_cores,
    target_cells,
    seed,
    iteration
) {
  selected_cells <- pca_stratified_cells(
    obj = obj,
    target_cells = target_cells,
    pcs_use = pcs_use,
    seed = seed
  )
  sweep_obj <- if (length(selected_cells) < ncol(obj)) {
    subset(obj, cells = selected_cells)
  } else {
    obj
  }
  sampling_method <- if (length(selected_cells) < ncol(obj)) {
    "pca_stratified"
  } else {
    "all_cells"
  }

  timed <- time_step(
    call_with_supported_args(
      fun = paramSweep_fun,
      args = list(
        seu = sweep_obj,
        PCs = seq_len(min(length(pcs_use), ncol(sweep_obj) - 1L)),
        sct = FALSE,
        pN = pN,
        num.cores = sweep_cores
      )
    )
  )
  sweep_stats <- summarizeSweep_fun(timed$value, GT = FALSE)
  selected <- select_best_pk(sweep_stats)
  evidence <- selected$bcmvn
  evidence$sweep_iteration <- iteration
  evidence$sweep_seed <- seed
  evidence$sweep_cells <- length(selected_cells)
  evidence$sampling_method <- sampling_method
  evidence$selected <- seq_len(nrow(evidence)) == selected$best_index

  rm(sweep_obj)
  gc()

  list(
    best_pk = selected$best_pk,
    best_bcmetric = selected$best_bcmetric,
    evidence = evidence,
    seconds = timed$seconds,
    sweep_cells = length(selected_cells),
    sampling_method = sampling_method
  )
}

resolve_pk <- function(
    obj,
    pcs_use,
    pN,
    sweep_cores,
    pk_mode,
    fixed_pk,
    sweep_min_cells,
    sweep_fraction,
    sweep_max_cells,
    sweep_seed_primary,
    sweep_seed_secondary,
    pk_tolerance,
    sweep_scale_factor
) {
  if (identical(pk_mode, "fixed")) {
    if (is.null(fixed_pk) || length(fixed_pk) != 1L ||
        is.na(fixed_pk) || !is.finite(fixed_pk) || fixed_pk <= 0) {
      stop("doublet.fixed_pk must be a positive finite number in fixed mode.")
    }

    return(list(
      best_pk = as.numeric(fixed_pk),
      best_bcmetric = NA_real_,
      evidence = data.frame(
        pK_numeric = as.numeric(fixed_pk),
        BCmetric_numeric = NA_real_,
        sweep_iteration = 0L,
        sweep_seed = NA_integer_,
        sweep_cells = 0L,
        sampling_method = "fixed",
        selected = TRUE
      ),
      seconds = 0,
      sweep_cells = 0L,
      sampling_method = "fixed",
      stable = TRUE,
      pk_delta = 0,
      attempts = 0L,
      source = "configured_fixed_pk"
    ))
  }

  if (identical(pk_mode, "sweep")) {
    result <- run_pk_sweep(
      obj, pcs_use, pN, sweep_cores, ncol(obj),
      sweep_seed_primary, 1L
    )
    return(c(
      result,
      list(
        stable = NA,
        pk_delta = NA_real_,
        attempts = 1L,
        source = "exhaustive_sweep"
      )
    ))
  }

  target <- resolve_sweep_size(
    ncol(obj), sweep_min_cells, sweep_fraction, sweep_max_cells
  )
  max_target <- min(as.integer(sweep_max_cells), ncol(obj))
  all_evidence <- list()
  total_seconds <- 0
  attempts <- 0L

  repeat {
    attempts <- attempts + 1L
    first <- run_pk_sweep(
      obj, pcs_use, pN, sweep_cores, target,
      sweep_seed_primary, attempts * 2L - 1L
    )
    second <- run_pk_sweep(
      obj, pcs_use, pN, sweep_cores, target,
      sweep_seed_secondary, attempts * 2L
    )
    all_evidence <- c(all_evidence, list(first$evidence, second$evidence))
    total_seconds <- total_seconds + first$seconds + second$seconds
    pk_delta <- abs(first$best_pk - second$best_pk)
    stable <- pk_delta <= pk_tolerance

    candidates <- list(first, second)
    metrics <- c(first$best_bcmetric, second$best_bcmetric)
    chosen <- candidates[[which.max(metrics)]]

    if (stable || target >= max_target) {
      evidence <- dplyr::bind_rows(all_evidence)
      evidence$selected_final <- FALSE
      final_rows <- which(
        evidence$selected &
          evidence$sweep_iteration %in% c(attempts * 2L - 1L, attempts * 2L) &
          abs(evidence$pK_numeric - chosen$best_pk) < .Machine$double.eps^0.5
      )
      if (length(final_rows) > 0) evidence$selected_final[[final_rows[[1]]]] <- TRUE

      return(list(
        best_pk = chosen$best_pk,
        best_bcmetric = chosen$best_bcmetric,
        evidence = evidence,
        seconds = total_seconds,
        sweep_cells = target,
        sampling_method = chosen$sampling_method,
        stable = stable,
        pk_delta = pk_delta,
        attempts = attempts * 2L,
        source = if (stable) "adaptive_stable" else "adaptive_max_cells"
      ))
    }

    expanded <- ceiling(target * sweep_scale_factor)
    target <- min(max_target, max(target + 1L, expanded))
  }
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
    status,
    pcs = NA_integer_,
    nfeatures = NA_integer_,
    sweep_cores = NA_integer_
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
    pk_mode = NA_character_,
    pk_source = status,
    pk_stable = NA,
    pk_delta = NA_real_,
    selected_bcmetric = NA_real_,
    sweep_cells = 0L,
    sweep_fraction = 0,
    sampling_method = "not_run",
    sweep_attempts = 0L,
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
    bcmvn = NULL,
    performance = data.frame(
      sample = sample_id,
      cells = length(cell_ids),
      nfeatures = nfeatures,
      pcs = pcs,
      sweep_cores = sweep_cores,
      pk_mode = NA_character_,
      pk_source = status,
      pk_stable = NA,
      pk_delta = NA_real_,
      selected_bcmetric = NA_real_,
      sweep_cells = 0L,
      sweep_fraction = 0,
      sampling_method = "not_run",
      sweep_attempts = 0L,
      normalize_seconds = NA_real_,
      find_variable_features_seconds = NA_real_,
      scale_data_seconds = NA_real_,
      run_pca_seconds = NA_real_,
      param_sweep_seconds = NA_real_,
      doublet_finder_seconds = NA_real_,
      total_seconds = 0,
      stringsAsFactors = FALSE
    )
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
    sweep_cores,
    pk_mode,
    fixed_pk,
    sweep_min_cells,
    sweep_fraction,
    sweep_max_cells,
    sweep_seed_primary,
    sweep_seed_secondary,
    pk_tolerance,
    sweep_scale_factor,
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
        status = "skipped_low_cell_number",
        pcs = pcs,
        nfeatures = nfeatures,
        sweep_cores = sweep_cores
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
        status = "skipped_too_few_pcs",
        pcs = max_pcs,
        nfeatures = nfeatures,
        sweep_cores = sweep_cores
      )
    )
  }

  message("[", sample_id, "] NormalizeData started")
  message("[", sample_id, "] FindVariableFeatures target: ", nfeatures)
  message("[", sample_id, "] ScaleData will use variable features only")
  message("[", sample_id, "] RunPCA target PCs: ", max_pcs)

  prep <- safe_scale_and_pca(
    obj = obj,
    npcs = max_pcs,
    nfeatures = nfeatures
  )
  obj <- prep$object
  variable_features <- prep$variable_features
  prep_timings <- prep$timings

  message_timed(
    sample_id,
    "NormalizeData",
    prep_timings$normalize_seconds
  )
  message_timed(
    sample_id,
    "FindVariableFeatures",
    prep_timings$find_variable_features_seconds
  )
  message_timed(
    sample_id,
    "ScaleData",
    prep_timings$scale_data_seconds
  )
  message_timed(
    sample_id,
    "RunPCA",
    prep_timings$run_pca_seconds
  )

  pcs_use <- seq_len(max_pcs)

  message(
    "[", sample_id, "] pK resolution started | mode = ", pk_mode,
    " | num.cores = ", sweep_cores
  )

  pk_res <- resolve_pk(
    obj = obj,
    pcs_use = pcs_use,
    pN = pN,
    sweep_cores = sweep_cores,
    pk_mode = pk_mode,
    fixed_pk = fixed_pk,
    sweep_min_cells = sweep_min_cells,
    sweep_fraction = sweep_fraction,
    sweep_max_cells = sweep_max_cells,
    sweep_seed_primary = sweep_seed_primary,
    sweep_seed_secondary = sweep_seed_secondary,
    pk_tolerance = pk_tolerance,
    sweep_scale_factor = sweep_scale_factor
  )
  best_pk <- pk_res$best_pk

  message_timed(sample_id, "pK resolution", pk_res$seconds)
  message(
    "[", sample_id, "] pK = ", best_pk,
    " | source = ", pk_res$source,
    " | stable = ", pk_res$stable,
    " | sweep cells = ", pk_res$sweep_cells,
    " / ", ncol(obj)
  )

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

  bcmvn <- pk_res$evidence
  bcmvn$sample <- sample_id
  bcmvn$pk_mode <- pk_mode

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

  message("[", sample_id, "] doubletFinder started")

  t_df <- time_step(
    call_with_supported_args(
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
  )
  obj <- t_df$value

  message_timed(
    sample_id,
    "doubletFinder",
    t_df$seconds
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
    pk_mode = pk_mode,
    pk_source = pk_res$source,
    pk_stable = pk_res$stable,
    pk_delta = pk_res$pk_delta,
    selected_bcmetric = pk_res$best_bcmetric,
    sweep_cells = pk_res$sweep_cells,
    sweep_fraction = pk_res$sweep_cells / ncol(obj),
    sampling_method = pk_res$sampling_method,
    sweep_attempts = pk_res$attempts,
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

  total_seconds <- sum(
    prep_timings$normalize_seconds,
    prep_timings$find_variable_features_seconds,
    prep_timings$scale_data_seconds,
    prep_timings$run_pca_seconds,
    pk_res$seconds,
    t_df$seconds,
    na.rm = TRUE
  )

  sample_performance <- data.frame(
    sample = sample_id,
    cells = ncol(obj),
    nfeatures = length(variable_features),
    pcs = length(pcs_use),
    sweep_cores = sweep_cores,
    pk_mode = pk_mode,
    pk_source = pk_res$source,
    pk_stable = pk_res$stable,
    pk_delta = pk_res$pk_delta,
    selected_bcmetric = pk_res$best_bcmetric,
    sweep_cells = pk_res$sweep_cells,
    sweep_fraction = pk_res$sweep_cells / ncol(obj),
    sampling_method = pk_res$sampling_method,
    sweep_attempts = pk_res$attempts,
    normalize_seconds = prep_timings$normalize_seconds,
    find_variable_features_seconds =
      prep_timings$find_variable_features_seconds,
    scale_data_seconds = prep_timings$scale_data_seconds,
    run_pca_seconds = prep_timings$run_pca_seconds,
    param_sweep_seconds = pk_res$seconds,
    doublet_finder_seconds = t_df$seconds,
    total_seconds = total_seconds,
    stringsAsFactors = FALSE
  )

  rm(obj)
  gc()

  list(
    metadata = sample_metadata,
    summary = sample_summary,
    bcmvn = bcmvn,
    performance = sample_performance
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
pk_mode <- tolower(as.character(
  get_config_value(config, "doublet.pk_mode", "auto")
))
fixed_pk <- get_config_value(
  config,
  "doublet.fixed_pk",
  get_config_value(config, "doublet.pK", NULL)
)
if (!is.null(fixed_pk)) fixed_pk <- suppressWarnings(as.numeric(fixed_pk))
sweep_min_cells <- as.integer(
  get_config_value(config, "doublet.sweep.min_cells", 2000)
)
sweep_fraction <- as.numeric(
  get_config_value(config, "doublet.sweep.fraction", 0.5)
)
sweep_max_cells <- as.integer(
  get_config_value(config, "doublet.sweep.max_cells", 5000)
)
sweep_seed_primary <- as.integer(
  get_config_value(config, "doublet.sweep.seed_primary", 1234)
)
sweep_seed_secondary <- as.integer(
  get_config_value(config, "doublet.sweep.seed_secondary", 4321)
)
pk_tolerance <- as.numeric(
  get_config_value(config, "doublet.sweep.pk_tolerance", 0.01)
)
sweep_scale_factor <- as.numeric(
  get_config_value(config, "doublet.sweep.scale_factor", 1.5)
)
sweep_cores <- if (!is.null(opt$cores)) {
  as.integer(opt$cores)
} else {
  as.integer(get_config_value(config, "doublet.sweep_cores", 1))
}

if (
  is.na(sweep_cores) ||
    length(sweep_cores) != 1 ||
    sweep_cores < 1
) {
  stop("doublet.sweep_cores must be a positive integer.")
}

if (!pk_mode %in% c("auto", "sweep", "fixed")) {
  stop("doublet.pk_mode must be one of: auto, sweep, fixed.")
}
if (is.na(sweep_min_cells) || sweep_min_cells < 1L) {
  stop("doublet.sweep.min_cells must be a positive integer.")
}
if (is.na(sweep_max_cells) || sweep_max_cells < sweep_min_cells) {
  stop("doublet.sweep.max_cells must be >= sweep.min_cells.")
}
if (is.na(sweep_fraction) || sweep_fraction <= 0 || sweep_fraction > 1) {
  stop("doublet.sweep.fraction must be in (0, 1].")
}
if (is.na(pk_tolerance) || pk_tolerance < 0) {
  stop("doublet.sweep.pk_tolerance must be >= 0.")
}
if (is.na(sweep_scale_factor) || sweep_scale_factor <= 1) {
  stop("doublet.sweep.scale_factor must be > 1.")
}

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
  sweep_cores = sweep_cores,
  pk_mode = pk_mode,
  fixed_pk = fixed_pk,
  sweep_min_cells = sweep_min_cells,
  sweep_fraction = sweep_fraction,
  sweep_max_cells = sweep_max_cells,
  sweep_seed_primary = sweep_seed_primary,
  sweep_seed_secondary = sweep_seed_secondary,
  pk_tolerance = pk_tolerance,
  sweep_scale_factor = sweep_scale_factor,
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

performance_table <- dplyr::bind_rows(
  lapply(results, `[[`, "performance")
)

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
    "sweep_cores",
    "pk_mode",
    "sweep_min_cells",
    "sweep_fraction",
    "sweep_max_cells",
    "sweep_seed_primary",
    "sweep_seed_secondary",
    "pk_tolerance",
    "sweep_scale_factor",
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
    sweep_cores,
    pk_mode,
    sweep_min_cells,
    sweep_fraction,
    sweep_max_cells,
    sweep_seed_primary,
    sweep_seed_secondary,
    pk_tolerance,
    sweep_scale_factor,
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

readr::write_csv(
  performance_table,
  "01g_doublet_detection_performance.csv"
)

sink("01g_sessionInfo.txt")
print(sessionInfo())
sink()

message("Doublet-scored object: ", opt$output)
message("No cells were removed by this module.")
message_step("01g_doublet_detection_scoring finished successfully")
