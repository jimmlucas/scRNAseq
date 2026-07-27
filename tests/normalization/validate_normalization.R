#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: validate_normalization.R ",
    "<input_doublet_filtered.rds> <normalization_output_dir>"
  )
}

input_object <- args[[1]]
output_dir <- args[[2]]

required_files <- c(
  "seu_normalized.rds",
  "04_normalization_summary.csv",
  "04_normalization_by_sample.csv",
  "04_variable_features.csv",
  "applied_normalization_policy.yaml",
  "04_normalization_sessionInfo.txt"
)

for (file_name in required_files) {
  path <- file.path(output_dir, file_name)
  if (!file.exists(path)) {
    stop("Missing expected output: ", path)
  }
}

if (!dir.exists(file.path(output_dir, "normalized_samples"))) {
  stop("Missing per-sample normalized object directory")
}

input <- readRDS(input_object)
normalized <- readRDS(file.path(output_dir, "seu_normalized.rds"))

if (!inherits(normalized, "Seurat")) {
  stop("Normalized output is not a Seurat object")
}

if (ncol(normalized) != ncol(input)) {
  stop("Cell count changed during normalization")
}

if (nrow(normalized) != nrow(input)) {
  stop("Feature count changed during normalization")
}

metadata <- normalized@meta.data

for (column_name in c(
  "normalization_method",
  "normalization_scope",
  "normalization_applied",
  "normalization_policy_status",
  "variable_features_method"
)) {
  if (!column_name %in% colnames(metadata)) {
    stop("Missing normalization metadata column: ", column_name)
  }
}

if (!all(metadata$normalization_applied %in% TRUE)) {
  stop("Not all cells are marked as normalized")
}

summary <- readr::read_csv(
  file.path(output_dir, "04_normalization_summary.csv"),
  show_col_types = FALSE
)
summary_map <- setNames(summary$value, summary$metric)

stopifnot(as.integer(summary_map[["cells_input"]]) == ncol(input))
stopifnot(as.integer(summary_map[["cells_output"]]) == ncol(input))
stopifnot(as.integer(summary_map[["features_input"]]) == nrow(input))
stopifnot(as.integer(summary_map[["features_output"]]) == nrow(input))
stopifnot(summary_map[["method"]] == "LogNormalize")
stopifnot(summary_map[["scope"]] == "per_sample")
stopifnot(summary_map[["assay"]] == "RNA")
stopifnot(summary_map[["application_status"]] == "applied")

by_sample <- readr::read_csv(
  file.path(output_dir, "04_normalization_by_sample.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(by_sample) == 2)
stopifnot(all(by_sample$method == "LogNormalize"))
stopifnot(all(by_sample$variable_features > 0))

variable_features <- readr::read_csv(
  file.path(output_dir, "04_variable_features.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(variable_features) > 0)
stopifnot(all(c("sample", "rank", "feature") %in% colnames(variable_features)))

sample_objects <- list.files(
  file.path(output_dir, "normalized_samples"),
  pattern = "\\.normalized\\.rds$",
  full.names = TRUE
)
stopifnot(length(sample_objects) == 2)

policy <- yaml::read_yaml(
  file.path(output_dir, "applied_normalization_policy.yaml")
)$normalization_policy

stopifnot(policy$approval$status == "approved")
stopifnot(policy$validation$status == "valid")
stopifnot(policy$application$status == "applied")
stopifnot(isTRUE(policy$application$normalization_applied))
stopifnot(as.integer(policy$application$cells_normalized) == ncol(input))
stopifnot(as.integer(policy$application$features_preserved) == nrow(input))

message("Normalization contract: VALID")
message("Cells: ", ncol(normalized))
message("Features: ", nrow(normalized))
message("Variable features rows: ", nrow(variable_features))

