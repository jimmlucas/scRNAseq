#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: validate_dimensionality.R ",
    "<input_normalized.rds> <dimensionality_output_dir>"
  )
}

input_object <- args[[1]]
output_dir <- args[[2]]

required_files <- c(
  "seu_pca.rds",
  "05_dimensionality_summary.csv",
  "05_variable_features_summary.csv",
  "05_pca_embeddings_preview.csv",
  "05_pca_loadings_preview.csv",
  "05_pca_elbow_metrics.csv",
  "05_dimensionality_sessionInfo.txt"
)

for (file_name in required_files) {
  path <- file.path(output_dir, file_name)
  if (!file.exists(path)) {
    stop("Missing expected output: ", path)
  }
}

input <- readRDS(input_object)
pca_object <- readRDS(file.path(output_dir, "seu_pca.rds"))

if (!inherits(pca_object, "Seurat")) {
  stop("PCA output is not a Seurat object")
}

if (ncol(pca_object) != ncol(input)) {
  stop("Cell count changed during dimensionality reduction")
}

if (nrow(pca_object) != nrow(input)) {
  stop("Feature count changed during dimensionality reduction")
}

if (!"pca" %in% names(pca_object@reductions)) {
  stop("Missing PCA reduction")
}

metadata <- pca_object@meta.data

for (column_name in c(
  "dimensionality_scaled",
  "dimensionality_reduction",
  "dimensionality_npcs",
  "dimensionality_assay",
  "dimensionality_policy_status"
)) {
  if (!column_name %in% colnames(metadata)) {
    stop("Missing dimensionality metadata column: ", column_name)
  }
}

if (!all(metadata$dimensionality_scaled %in% TRUE)) {
  stop("Not all cells are marked as scaled")
}

if (!all(metadata$dimensionality_reduction == "pca")) {
  stop("Unexpected dimensionality reduction metadata")
}

summary <- readr::read_csv(
  file.path(output_dir, "05_dimensionality_summary.csv"),
  show_col_types = FALSE
)
summary_map <- setNames(summary$value, summary$metric)

stopifnot(as.integer(summary_map[["cells_input"]]) == ncol(input))
stopifnot(as.integer(summary_map[["cells_output"]]) == ncol(input))
stopifnot(as.integer(summary_map[["features_input"]]) == nrow(input))
stopifnot(as.integer(summary_map[["features_output"]]) == nrow(input))
stopifnot(summary_map[["scaling_applied"]] == "true")
stopifnot(summary_map[["pca_applied"]] == "true")
stopifnot(summary_map[["normalization_application_status"]] == "applied")

npcs <- as.integer(summary_map[["npcs_computed"]])
stopifnot(npcs >= 2)
stopifnot(ncol(Embeddings(pca_object, "pca")) == npcs)

variable_features <- readr::read_csv(
  file.path(output_dir, "05_variable_features_summary.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(variable_features) > 0)
stopifnot(all(c("rank", "feature") %in% colnames(variable_features)))

embeddings_preview <- readr::read_csv(
  file.path(output_dir, "05_pca_embeddings_preview.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(embeddings_preview) > 0)
stopifnot("cell" %in% colnames(embeddings_preview))
stopifnot("PC_1" %in% colnames(embeddings_preview))

loadings_preview <- readr::read_csv(
  file.path(output_dir, "05_pca_loadings_preview.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(loadings_preview) > 0)
stopifnot("feature" %in% colnames(loadings_preview))
stopifnot("PC_1" %in% colnames(loadings_preview))

elbow <- readr::read_csv(
  file.path(output_dir, "05_pca_elbow_metrics.csv"),
  show_col_types = FALSE
)
stopifnot(nrow(elbow) == npcs)
stopifnot(all(c(
  "pc",
  "stdev",
  "variance",
  "variance_percent",
  "cumulative_variance_percent"
) %in% colnames(elbow)))

message("Dimensionality contract: VALID")
message("Cells: ", ncol(pca_object))
message("Features: ", nrow(pca_object))
message("PCs: ", npcs)
