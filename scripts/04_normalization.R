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
    "--policy",
    type = "character",
    default = NULL,
    help = "Validated and approved normalization policy YAML"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_normalized.rds",
    help = "Output normalized Seurat RDS object"
  ),
  make_option(
    "--sample-output-dir",
    dest = "sample_output_dir",
    type = "character",
    default = "normalized_samples",
    help = "Directory for per-sample normalized Seurat objects"
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

require_positive_number <- function(value, label) {
  if (is.null(value) || length(value) != 1 || is.na(value)) {
    stop(label, " must be a non-missing positive number")
  }

  parsed <- suppressWarnings(as.numeric(value))

  if (is.na(parsed) || !is.finite(parsed) || parsed <= 0) {
    stop(label, " must be a positive finite number")
  }

  parsed
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

safe_file_token <- function(value) {
  token <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  token <- gsub("^_+|_+$", "", token)

  if (!nzchar(token)) {
    return("sample")
  }

  token
}

normalise_regressors <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(character())
  }

  if (length(value) == 1) {
    text <- tolower(as.character(value))
    if (!nzchar(text) || text %in% c("none", "null")) {
      return(character())
    }
  }

  unique(as.character(value))
}

check_file(opt$input, "Doublet-filtered Seurat object")
check_file(opt$policy, "Approved normalization policy")

message_step("Loading approved normalization policy")

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

if (!isTRUE(validation$validated_for_normalization)) {
  stop(
    "Normalization policy has not been validated for NORMALIZATION"
  )
}

if (!identical(application_status, "not_applied")) {
  stop(
    "Normalization policy application.status must be not_applied. Observed: ",
    application_status
  )
}

if (!identical(require_bool(
  application$normalization_applied,
  "application.normalization_applied"
), FALSE)) {
  stop("application.normalization_applied must be false")
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
scope <- require_scalar(strategy$scope, "strategy.scope")
assay <- require_scalar(strategy$assay, "strategy.assay")
variable_features_method <- require_scalar(
  strategy$variable_features_method,
  "strategy.variable_features_method"
)
nfeatures <- require_positive_int(
  strategy$nfeatures,
  "strategy.nfeatures"
)

if (!method %in% c("LogNormalize", "SCTransform")) {
  stop("Unsupported normalization method: ", method)
}

if (!identical(scope, "per_sample")) {
  stop("Only per_sample normalization scope is supported in this MVP")
}

if (!identical(assay, "RNA")) {
  stop("Only RNA assay normalization is supported in this MVP")
}

if (!identical(variable_features_method, "vst")) {
  stop("Only vst variable feature selection is supported in this MVP")
}

scale_factor <- NA_real_
regressors <- normalise_regressors(strategy$regressors)

if (identical(method, "LogNormalize")) {
  scale_factor <- require_positive_number(
    strategy$scale_factor,
    "strategy.scale_factor"
  )

  if (length(regressors) > 0) {
    stop("LogNormalize policy must keep strategy.regressors as none")
  }
}

message("Project ID: ", policy_project_id)
message("Biological run ID: ", policy_run_id)
message("Method: ", method)
message("Scope: ", scope)

message_step("Loading Seurat object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
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
expected_samples <- require_positive_int(input_contract$samples, "input.samples")

if (cells_input != expected_cells) {
  stop(
    "Input cell count does not match policy: ",
    cells_input,
    " != ",
    expected_cells
  )
}

if (features_input != expected_features) {
  stop(
    "Input feature count does not match policy: ",
    features_input,
    " != ",
    expected_features
  )
}

sample_column <- find_sample_column(seu@meta.data)
sample_values <- as.character(seu@meta.data[[sample_column]])
sample_names <- sort(unique(sample_values))

if (length(sample_names) != expected_samples) {
  stop(
    "Input sample count does not match policy: ",
    length(sample_names),
    " != ",
    expected_samples
  )
}

dir.create(opt$sample_output_dir, recursive = TRUE, showWarnings = FALSE)

message_step("Normalizing samples")

sample_objects <- list()
sample_summary <- list()
variable_feature_rows <- list()

for (sample_name in sample_names) {
  sample_cells <- rownames(seu@meta.data)[sample_values == sample_name]

  message(
    "Normalizing sample ",
    sample_name,
    " (",
    length(sample_cells),
    " cells)"
  )

  sample_obj <- subset(seu, cells = sample_cells)
  DefaultAssay(sample_obj) <- assay

  if (identical(method, "LogNormalize")) {
    sample_obj <- NormalizeData(
      sample_obj,
      assay = assay,
      normalization.method = "LogNormalize",
      scale.factor = scale_factor,
      verbose = FALSE
    )

    sample_obj <- FindVariableFeatures(
      sample_obj,
      assay = assay,
      selection.method = variable_features_method,
      nfeatures = nfeatures,
      verbose = FALSE
    )
  } else if (identical(method, "SCTransform")) {
    vars_to_regress <- regressors

    sample_obj <- SCTransform(
      sample_obj,
      assay = assay,
      variable.features.n = nfeatures,
      vars.to.regress = vars_to_regress,
      verbose = FALSE
    )
  }

  sample_obj$normalization_method <- method
  sample_obj$normalization_scope <- scope
  sample_obj$normalization_assay <- assay
  sample_obj$normalization_applied <- TRUE
  sample_obj$normalization_policy_status <- approval_status
  sample_obj$variable_features_method <- variable_features_method

  variable_features <- VariableFeatures(sample_obj)
  variable_features <- head(variable_features, nfeatures)

  sample_token <- safe_file_token(sample_name)
  sample_path <- file.path(
    opt$sample_output_dir,
    paste0(sample_token, ".normalized.rds")
  )
  saveRDS(sample_obj, sample_path)

  sample_objects[[sample_name]] <- sample_obj

  sample_summary[[sample_name]] <- tibble::tibble(
    sample = sample_name,
    cells = ncol(sample_obj),
    features = nrow(sample_obj),
    method = method,
    assay = DefaultAssay(sample_obj),
    variable_features_method = variable_features_method,
    variable_features = length(variable_features),
    nfeatures_requested = nfeatures,
    scale_factor = ifelse(is.na(scale_factor), "", as.character(scale_factor)),
    regressors = ifelse(
      length(regressors) == 0,
      "none",
      paste(regressors, collapse = ",")
    ),
    output_file = sample_path
  )

  if (length(variable_features) > 0) {
    variable_feature_rows[[sample_name]] <- tibble::tibble(
      sample = sample_name,
      rank = seq_along(variable_features),
      feature = variable_features
    )
  } else {
    variable_feature_rows[[sample_name]] <- tibble::tibble(
      sample = character(),
      rank = integer(),
      feature = character()
    )
  }

  rm(sample_obj)
  gc(verbose = FALSE)
}

message_step("Combining normalized samples")

combined <- sample_objects[[1]]

if (length(sample_objects) > 1) {
  combined <- merge(
    sample_objects[[1]],
    y = sample_objects[-1],
    merge.data = TRUE
  )
}

DefaultAssay(combined) <- ifelse(
  identical(method, "SCTransform") && "SCT" %in% names(combined@assays),
  "SCT",
  assay
)

combined$normalization_method <- method
combined$normalization_scope <- scope
combined$normalization_assay <- assay
combined$normalization_applied <- TRUE
combined$normalization_policy_status <- approval_status
combined$variable_features_method <- variable_features_method

variable_features_table <- dplyr::bind_rows(variable_feature_rows)

if (nrow(variable_features_table) > 0) {
  combined_variable_features <- unique(variable_features_table$feature)
  combined_variable_features <- combined_variable_features[
    combined_variable_features %in% rownames(combined)
  ]
  VariableFeatures(combined) <- combined_variable_features
}

cells_output <- ncol(combined)
features_output <- nrow(combined)

if (cells_output != cells_input) {
  stop(
    "Normalization changed cell count: ",
    cells_output,
    " != ",
    cells_input
  )
}

if (features_output != features_input) {
  stop(
    "Normalization changed feature count: ",
    features_output,
    " != ",
    features_input
  )
}

saveRDS(combined, opt$output)

message_step("Writing audit artifacts")

summary_table <- tibble::tibble(
  metric = c(
    "cells_input",
    "cells_output",
    "features_input",
    "features_output",
    "samples",
    "method",
    "scope",
    "assay",
    "variable_features_method",
    "nfeatures_requested",
    "scale_factor",
    "regressors",
    "approval_status",
    "validation_status",
    "application_status",
    "normalized_object"
  ),
  value = c(
    as.character(cells_input),
    as.character(cells_output),
    as.character(features_input),
    as.character(features_output),
    as.character(length(sample_names)),
    method,
    scope,
    assay,
    variable_features_method,
    as.character(nfeatures),
    ifelse(is.na(scale_factor), "", as.character(scale_factor)),
    ifelse(length(regressors) == 0, "none", paste(regressors, collapse = ",")),
    approval_status,
    validation_status,
    "applied",
    opt$output
  )
)

readr::write_csv(
  summary_table,
  "04_normalization_summary.csv"
)

readr::write_csv(
  dplyr::bind_rows(sample_summary),
  "04_normalization_by_sample.csv"
)

readr::write_csv(
  variable_features_table,
  "04_variable_features.csv"
)

applied_at <- format(
  Sys.time(),
  "%Y-%m-%dT%H:%M:%S%z"
)
applied_at <- sub("(\\+|\\-)([0-9]{2})([0-9]{2})$", "\\1\\2:\\3", applied_at)

applied_policy <- policy_document
applied_policy$normalization_policy$application$status <- "applied"
applied_policy$normalization_policy$application$normalization_applied <- TRUE
applied_policy$normalization_policy$application$cells_normalized <- cells_output
applied_policy$normalization_policy$application$features_preserved <- features_output
applied_policy$normalization_policy$application$samples_normalized <- length(sample_names)
applied_policy$normalization_policy$application$applied_at <- applied_at
applied_policy$normalization_policy$application$output_object <- opt$output
applied_policy$normalization_policy$application$sample_output_dir <-
  opt$sample_output_dir

yaml::write_yaml(
  applied_policy,
  "applied_normalization_policy.yaml"
)

sink("04_normalization_sessionInfo.txt")
print(sessionInfo())
sink()

message("Normalized object saved in: ", opt$output)
message("Per-sample objects saved in: ", opt$sample_output_dir)
message("Summary saved in: 04_normalization_summary.csv")
