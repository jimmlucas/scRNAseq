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
    help = "Normalized Seurat RDS object"
  ),
  make_option(
    "--policy",
    type = "character",
    default = NULL,
    help = "Applied normalization policy YAML"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_pca.rds",
    help = "Output PCA-ready Seurat RDS object"
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
  ),
  make_option(
    "--npcs",
    type = "integer",
    default = 30,
    help = "Requested number of principal components"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

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

require_bool <- function(value, label) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(label, " must be a non-missing boolean")
  }

  value
}

require_positive_int <- function(value, label) {
  if (is.null(value) || length(value) != 1 || is.na(value)) {
    stop(label, " must be a non-missing positive integer")
  }

  parsed <- suppressWarnings(as.numeric(value))

  if (is.na(parsed) || parsed <= 0 || parsed != floor(parsed)) {
    stop(label, " must be a positive integer")
  }

  as.integer(parsed)
}

check_file(opt$input, "Normalized Seurat object")
check_file(opt$policy, "Applied normalization policy")

if (is.null(opt$npcs) || is.na(opt$npcs) || opt$npcs <= 0) {
  stop("--npcs must be a positive integer")
}

message_step("Loading applied normalization policy")

policy_document <- yaml::read_yaml(opt$policy)
policy <- require_mapping(
  policy_document$normalization_policy,
  "normalization_policy"
)

identity <- require_mapping(
  policy$identity,
  "normalization_policy.identity"
)
input_contract <- require_mapping(
  policy$input,
  "normalization_policy.input"
)
strategy <- require_mapping(
  policy$strategy,
  "normalization_policy.strategy"
)
approval <- require_mapping(
  policy$approval,
  "normalization_policy.approval"
)
validation <- require_mapping(
  policy$validation,
  "normalization_policy.validation"
)
application <- require_mapping(
  policy$application,
  "normalization_policy.application"
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

if (!approval_status %in% c("approved", "modified_and_approved")) {
  stop(
    "Normalization policy approval.status must be approved or ",
    "modified_and_approved. Observed: ",
    approval_status
  )
}

if (!identical(validation_status, "valid")) {
  stop(
    "Normalization policy validation.status must be valid. Observed: ",
    validation_status
  )
}

if (!identical(application_status, "applied")) {
  stop(
    "Normalization policy application.status must be applied. Observed: ",
    application_status
  )
}

if (!identical(require_bool(
  application$normalization_applied,
  "application.normalization_applied"
), TRUE)) {
  stop("application.normalization_applied must be true")
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
  stop(
    "Policy project_id does not match execution metadata: ",
    policy_project_id,
    " != ",
    opt$project_id
  )
}

if (
  nzchar(opt$biological_run_id) &&
  policy_run_id != opt$biological_run_id
) {
  stop(
    "Policy biological_run_id does not match execution metadata: ",
    policy_run_id,
    " != ",
    opt$biological_run_id
  )
}

method <- require_scalar(strategy$method, "strategy.method")
assay <- require_scalar(strategy$assay, "strategy.assay")
nfeatures <- require_positive_int(
  strategy$nfeatures,
  "strategy.nfeatures"
)

message("Project ID: ", policy_project_id)
message("Biological run ID: ", policy_run_id)
message("Normalization method: ", method)

message_step("Loading normalized Seurat object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

if (identical(method, "SCTransform") && "SCT" %in% names(seu@assays)) {
  assay <- "SCT"
}

if (!assay %in% names(seu@assays)) {
  stop("Assay not found in Seurat object: ", assay)
}

DefaultAssay(seu) <- assay

cells_input <- ncol(seu)
features_input <- nrow(seu)

expected_cells <- require_positive_int(input_contract$cells, "input.cells")
expected_features <- require_positive_int(
  input_contract$features,
  "input.features"
)

if (cells_input != expected_cells) {
  stop(
    "Input cell count does not match normalization policy: ",
    cells_input,
    " != ",
    expected_cells
  )
}

if (features_input != expected_features) {
  stop(
    "Input feature count does not match normalization policy: ",
    features_input,
    " != ",
    expected_features
  )
}

variable_features <- VariableFeatures(seu)
variable_features <- variable_features[
  variable_features %in% rownames(seu)
]

if (length(variable_features) == 0) {
  stop("Input object has no variable features for PCA")
}

if (length(variable_features) < 2) {
  stop("At least two variable features are required for PCA")
}

max_pcs <- min(
  as.integer(opt$npcs),
  length(variable_features),
  cells_input - 1
)

if (max_pcs < 2) {
  stop(
    "At least two principal components are required. Computable PCs: ",
    max_pcs
  )
}

message_step("Scaling variable features")

seu <- ScaleData(
  seu,
  assay = assay,
  features = variable_features,
  verbose = FALSE
)

message_step("Running PCA")

seu <- RunPCA(
  seu,
  assay = assay,
  features = variable_features,
  npcs = max_pcs,
  verbose = FALSE
)

seu$dimensionality_scaled <- TRUE
seu$dimensionality_reduction <- "pca"
seu$dimensionality_npcs <- max_pcs
seu$dimensionality_assay <- assay
seu$dimensionality_policy_status <- application_status

cells_output <- ncol(seu)
features_output <- nrow(seu)

if (cells_output != cells_input) {
  stop("PCA changed cell count: ", cells_output, " != ", cells_input)
}

if (features_output != features_input) {
  stop("PCA changed feature count: ", features_output, " != ", features_input)
}

saveRDS(seu, opt$output)

message_step("Writing audit artifacts")

pca <- seu[["pca"]]
stdev <- pca@stdev
variance <- stdev^2
variance_percent <- variance / sum(variance) * 100
cumulative_variance_percent <- cumsum(variance_percent)
pc_names <- paste0("PC_", seq_along(stdev))

summary_table <- tibble::tibble(
  metric = c(
    "cells_input",
    "cells_output",
    "features_input",
    "features_output",
    "assay",
    "normalization_method",
    "variable_features",
    "nfeatures_requested",
    "npcs_requested",
    "npcs_computed",
    "scaling_applied",
    "pca_applied",
    "approval_status",
    "validation_status",
    "normalization_application_status",
    "pca_object"
  ),
  value = c(
    as.character(cells_input),
    as.character(cells_output),
    as.character(features_input),
    as.character(features_output),
    assay,
    method,
    as.character(length(variable_features)),
    as.character(nfeatures),
    as.character(opt$npcs),
    as.character(max_pcs),
    "true",
    "true",
    approval_status,
    validation_status,
    application_status,
    opt$output
  )
)

readr::write_csv(
  summary_table,
  "05_dimensionality_summary.csv"
)

readr::write_csv(
  tibble::tibble(
    rank = seq_along(variable_features),
    feature = variable_features
  ),
  "05_variable_features_summary.csv"
)

embeddings <- as.data.frame(Embeddings(seu, "pca"))
embeddings$cell <- rownames(embeddings)
preview_pcs <- head(colnames(embeddings)[colnames(embeddings) != "cell"], 10)
embedding_preview <- embeddings[
  seq_len(min(nrow(embeddings), 100)),
  c("cell", preview_pcs),
  drop = FALSE
]

readr::write_csv(
  embedding_preview,
  "05_pca_embeddings_preview.csv"
)

loadings <- as.data.frame(Loadings(seu, "pca"))
loadings$feature <- rownames(loadings)
loading_preview_pcs <- head(
  colnames(loadings)[colnames(loadings) != "feature"],
  10
)
loading_preview <- loadings[
  seq_len(min(nrow(loadings), 100)),
  c("feature", loading_preview_pcs),
  drop = FALSE
]

readr::write_csv(
  loading_preview,
  "05_pca_loadings_preview.csv"
)

readr::write_csv(
  tibble::tibble(
    pc = pc_names,
    stdev = stdev,
    variance = variance,
    variance_percent = variance_percent,
    cumulative_variance_percent = cumulative_variance_percent
  ),
  "05_pca_elbow_metrics.csv"
)

sink("05_dimensionality_sessionInfo.txt")
print(sessionInfo())
sink()

message("Dimensionality reduction completed successfully")
