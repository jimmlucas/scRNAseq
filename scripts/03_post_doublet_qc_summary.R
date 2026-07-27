#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

option_list <- list(
  make_option(
    "--input",
    type = "character",
    default = NULL,
    help = "Doublet-filtered Seurat RDS object"
  ),
  make_option(
    "--applied-policy",
    dest = "applied_policy",
    type = "character",
    default = NULL,
    help = "Applied doublet policy YAML"
  ),
  make_option(
    "--min-cells-per-sample",
    dest = "min_cells_per_sample",
    type = "integer",
    default = 50,
    help = "Warning threshold for final cells per sample"
  ),
  make_option(
    "--project-id",
    dest = "project_id",
    type = "character",
    default = "",
    help = "Expected stable project identifier"
  ),
  make_option(
    "--biological-run-id",
    dest = "biological_run_id",
    type = "character",
    default = "",
    help = "Expected stable biological run identifier"
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

require_mapping <- function(value, label) {
  if (is.null(value) || !is.list(value)) {
    stop(label, " must be a YAML mapping")
  }

  value
}

require_scalar <- function(value, label) {
  if (is.null(value) || length(value) != 1 || is.na(value)) {
    stop(label, " must be a non-missing scalar")
  }

  as.character(value)
}

require_int <- function(value, label) {
  if (is.null(value) || length(value) != 1 || is.na(value)) {
    stop(label, " must be a non-missing integer")
  }

  parsed <- suppressWarnings(as.integer(value))

  if (is.na(parsed) || parsed < 0) {
    stop(label, " must be a non-negative integer")
  }

  parsed
}

safe_numeric_summary <- function(values, prefix) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]

  if (length(finite_values) == 0) {
    return(
      tibble::tibble(
        metric = c(
          paste0(prefix, "_min"),
          paste0(prefix, "_q25"),
          paste0(prefix, "_median"),
          paste0(prefix, "_mean"),
          paste0(prefix, "_q75"),
          paste0(prefix, "_max")
        ),
        value = NA_real_
      )
    )
  }

  tibble::tibble(
    metric = c(
      paste0(prefix, "_min"),
      paste0(prefix, "_q25"),
      paste0(prefix, "_median"),
      paste0(prefix, "_mean"),
      paste0(prefix, "_q75"),
      paste0(prefix, "_max")
    ),
    value = c(
      min(finite_values),
      unname(stats::quantile(finite_values, 0.25)),
      stats::median(finite_values),
      mean(finite_values),
      unname(stats::quantile(finite_values, 0.75)),
      max(finite_values)
    )
  )
}

find_sample_column <- function(metadata) {
  for (column_name in c("sample", "ShortName", "sample_id", "orig.ident")) {
    if (column_name %in% colnames(metadata)) {
      return(column_name)
    }
  }

  stop(
    paste(
      "No sample metadata column found.",
      "Expected one of: sample, ShortName, sample_id, orig.ident"
    )
  )
}

detect_metric_column <- function(metadata, candidates) {
  found <- candidates[candidates %in% colnames(metadata)]

  if (length(found) == 0) {
    return(NULL)
  }

  found[[1]]
}

append_warning <- function(
    warnings,
    severity,
    code,
    message,
    sample = "",
    metric = "",
    value = NA,
    threshold = NA
) {
  dplyr::bind_rows(
    warnings,
    tibble::tibble(
      severity = severity,
      code = code,
      message = message,
      sample = sample,
      metric = metric,
      value = as.character(value),
      threshold = as.character(threshold)
    )
  )
}

check_file(opt$input, "Doublet-filtered Seurat object")
check_file(opt$applied_policy, "Applied doublet policy")

message_step("Loading applied doublet policy")

policy_document <- yaml::read_yaml(opt$applied_policy)
policy <- require_mapping(
  policy_document$doublet_policy,
  "doublet_policy"
)
identity <- require_mapping(
  policy$identity,
  "doublet_policy.identity"
)
approval <- require_mapping(
  policy$approval,
  "doublet_policy.approval"
)
validation <- require_mapping(
  policy$validation,
  "doublet_policy.validation"
)
application <- require_mapping(
  policy$application,
  "doublet_policy.application"
)
sample_policies <- require_mapping(
  policy$samples,
  "doublet_policy.samples"
)

approval_status <- tolower(
  require_scalar(approval$status, "approval.status")
)
validation_status <- tolower(
  require_scalar(validation$status, "validation.status")
)
application_status <- tolower(
  require_scalar(application$status, "application.status")
)

if (!identical(approval_status, "approved")) {
  stop("approval.status must be approved")
}

if (!identical(validation_status, "valid")) {
  stop("validation.status must be valid")
}

if (!identical(application_status, "applied")) {
  stop("application.status must be applied")
}

if (!isTRUE(application$filtering_applied)) {
  stop("application.filtering_applied must be true")
}

policy_project_id <- require_scalar(
  identity$project_id,
  "identity.project_id"
)
policy_run_id <- require_scalar(
  identity$biological_run_id,
  "identity.biological_run_id"
)

if (nzchar(opt$project_id) && policy_project_id != opt$project_id) {
  stop("Policy project_id does not match execution metadata")
}

if (
  nzchar(opt$biological_run_id) &&
  policy_run_id != opt$biological_run_id
) {
  stop("Policy biological_run_id does not match execution metadata")
}

message_step("Loading post-doublet object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

cells_final <- ncol(seu)
features_final <- nrow(seu)

if (cells_final == 0) {
  stop("Post-doublet object contains no cells")
}

metadata <- seu@meta.data

if (!identical(rownames(metadata), colnames(seu))) {
  stop("Metadata row names do not match Seurat cell IDs")
}

sample_column <- find_sample_column(metadata)
sample_values <- as.character(metadata[[sample_column]])

policy_samples <- sort(names(sample_policies))
object_samples <- sort(unique(sample_values))

if (!setequal(policy_samples, object_samples)) {
  stop(
    "Applied policy samples do not match filtered object samples. policy=",
    paste(policy_samples, collapse = ","),
    "; object=",
    paste(object_samples, collapse = ",")
  )
}

required_filter_columns <- c(
  "doublet_filter_action",
  "doublet_filter_applied",
  "doublet_filter_removed",
  "doublet_filter_policy_status"
)

missing_filter_columns <- setdiff(
  required_filter_columns,
  colnames(metadata)
)

if (length(missing_filter_columns) > 0) {
  stop(
    "Missing post-doublet filter metadata columns: ",
    paste(missing_filter_columns, collapse = ", ")
  )
}

if (any(metadata$doublet_filter_removed, na.rm = TRUE)) {
  stop(
    "Filtered object still contains cells marked as removed by DOUBLET_FILTER"
  )
}

metric_nfeature <- detect_metric_column(
  metadata,
  c("nFeature_RNA")
)
metric_ncount <- detect_metric_column(
  metadata,
  c("nCount_RNA")
)
metric_mito <- detect_metric_column(
  metadata,
  c("percent.mito", "percent.mt", "percent_mito", "percent_mt")
)

warnings <- tibble::tibble(
  severity = character(),
  code = character(),
  message = character(),
  sample = character(),
  metric = character(),
  value = character(),
  threshold = character()
)

if (is.null(metric_nfeature)) {
  warnings <- append_warning(
    warnings,
    "warning",
    "MISSING_NFEATURE_RNA",
    "nFeature_RNA is not available for post-doublet QC summary.",
    metric = "nFeature_RNA"
  )
}

if (is.null(metric_ncount)) {
  warnings <- append_warning(
    warnings,
    "warning",
    "MISSING_NCOUNT_RNA",
    "nCount_RNA is not available for post-doublet QC summary.",
    metric = "nCount_RNA"
  )
}

if (is.null(metric_mito)) {
  warnings <- append_warning(
    warnings,
    "warning",
    "MISSING_MITO_PERCENT",
    "Mitochondrial percentage metric is not available.",
    metric = "percent_mito"
  )
}

message_step("Checking per-sample status")

by_sample_rows <- lapply(
  policy_samples,
  function(sample_name) {
    sample_policy <- require_mapping(
      sample_policies[[sample_name]],
      paste0("samples.", sample_name)
    )

    sample_index <- which(sample_values == sample_name)
    final_cells <- length(sample_index)
    original_cells <- require_int(
      sample_policy$cells_scored,
      paste0("samples.", sample_name, ".cells_scored")
    )
    approved_doublets <- require_int(
      sample_policy$proposed_doublet_count,
      paste0("samples.", sample_name, ".proposed_doublet_count")
    )
    action <- require_scalar(
      sample_policy$application_action,
      paste0("samples.", sample_name, ".application_action")
    )
    expected_removed <- if (identical(action, "remove_proposed_doublets")) {
      approved_doublets
    } else {
      0L
    }
    expected_final <- original_cells - expected_removed

    if (final_cells != expected_final) {
      stop(
        "Final cell count mismatch for sample ",
        sample_name,
        ": expected=",
        expected_final,
        ", observed=",
        final_cells
      )
    }

    if (final_cells < opt$min_cells_per_sample) {
      warnings <<- append_warning(
        warnings,
        "warning",
        "LOW_FINAL_CELL_COUNT",
        "Sample has fewer cells than the configured readiness threshold.",
        sample = sample_name,
        metric = "cells_final",
        value = final_cells,
        threshold = opt$min_cells_per_sample
      )
    }

    sample_metadata <- metadata[sample_index, , drop = FALSE]

    tibble::tibble(
      sample = sample_name,
      action = action,
      cells_original = original_cells,
      doublets_removed = expected_removed,
      cells_final = final_cells,
      median_nFeature_RNA = if (!is.null(metric_nfeature)) {
        stats::median(
          as.numeric(sample_metadata[[metric_nfeature]]),
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      median_nCount_RNA = if (!is.null(metric_ncount)) {
        stats::median(
          as.numeric(sample_metadata[[metric_ncount]]),
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      median_percent_mito = if (!is.null(metric_mito)) {
        stats::median(
          as.numeric(sample_metadata[[metric_mito]]),
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
    )
  }
)

by_sample <- dplyr::bind_rows(by_sample_rows)

cells_from_policy <- sum(
  vapply(
    sample_policies,
    function(sample_policy) {
      require_int(sample_policy$cells_scored, "samples.cells_scored")
    },
    integer(1)
  )
)
doublets_removed <- require_int(
  application$cells_removed,
  "application.cells_removed"
)
expected_cells_final <- cells_from_policy - doublets_removed

if (cells_final != expected_cells_final) {
  stop(
    "Final cell count does not match applied policy: expected=",
    expected_cells_final,
    ", observed=",
    cells_final
  )
}

if (!is.null(metric_nfeature)) {
  if (anyNA(metadata[[metric_nfeature]])) {
    warnings <- append_warning(
      warnings,
      "warning",
      "MISSING_NFEATURE_VALUES",
      "nFeature_RNA contains missing values.",
      metric = "nFeature_RNA"
    )
  }
}

if (!is.null(metric_ncount)) {
  if (anyNA(metadata[[metric_ncount]])) {
    warnings <- append_warning(
      warnings,
      "warning",
      "MISSING_NCOUNT_VALUES",
      "nCount_RNA contains missing values.",
      metric = "nCount_RNA"
    )
  }
}

assays_available <- names(seu@assays)
modalities_available <- assays_available

message_step("Creating readiness outputs")

summary_rows <- dplyr::bind_rows(
  tibble::tibble(
    metric = c(
      "cells",
      "features",
      "samples",
      "doublets_removed",
      "min_cells_per_sample",
      "warning_count"
    ),
    value = c(
      cells_final,
      features_final,
      length(policy_samples),
      doublets_removed,
      opt$min_cells_per_sample,
      nrow(warnings)
    )
  ),
  if (!is.null(metric_nfeature)) {
    safe_numeric_summary(
      metadata[[metric_nfeature]],
      "nFeature_RNA"
    )
  },
  if (!is.null(metric_ncount)) {
    safe_numeric_summary(
      metadata[[metric_ncount]],
      "nCount_RNA"
    )
  },
  if (!is.null(metric_mito)) {
    safe_numeric_summary(
      metadata[[metric_mito]],
      "percent_mito"
    )
  }
)

readiness_status <- if (any(warnings$severity == "error")) {
  "blocked"
} else {
  "ready"
}

readiness <- list(
  analysis_readiness = list(
    status = readiness_status,
    qc = list(
      completed = TRUE
    ),
    doublets = list(
      scoring_completed = TRUE,
      human_approval_validated = TRUE,
      filtering_applied = TRUE
    ),
    dataset = list(
      cells = cells_final,
      features = features_final,
      samples = length(policy_samples),
      modalities = modalities_available
    ),
    warnings = list(
      count = nrow(warnings)
    ),
    next_phase = list(
      recommended = "normalization_policy_resolution"
    )
  )
)

message_step("Saving outputs")

readr::write_csv(
  summary_rows,
  "03_post_doublet_qc_summary.csv"
)

readr::write_csv(
  by_sample,
  "03_post_doublet_qc_by_sample.csv"
)

readr::write_tsv(
  warnings,
  "03_post_doublet_qc_warnings.tsv"
)

yaml::write_yaml(
  readiness,
  "03_analysis_readiness.yaml"
)

sink("03_sessionInfo.txt")
print(sessionInfo())
sink()

message("Analysis readiness: ", readiness_status)
message("Cells: ", cells_final)
message("Features: ", features_final)
message("Samples: ", length(policy_samples))
