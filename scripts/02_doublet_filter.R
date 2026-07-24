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
    help = "Doublet-scored Seurat RDS object"
  ),
  make_option(
    "--policy",
    type = "character",
    default = NULL,
    help = "Validated and approved doublet policy YAML"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_doublet_filtered.rds",
    help = "Output doublet-filtered Seurat RDS object"
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

allowed_actions <- c(
  "remove_proposed_doublets",
  "keep_all",
  "exclude_sample_from_filtering"
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

normalise_logical <- function(value) {
  if (is.logical(value)) {
    return(value)
  }

  if (is.numeric(value)) {
    return(value == 1)
  }

  tolower(as.character(value)) %in% c("true", "t", "1", "yes")
}

check_file(opt$input, "Input Seurat object")
check_file(opt$policy, "Approved doublet policy")

message_step("Loading approved doublet policy")

policy_document <- yaml::read_yaml(opt$policy)
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
  stop(
    "Doublet policy approval.status must be approved. Observed: ",
    approval_status
  )
}

if (!identical(validation_status, "valid")) {
  stop(
    "Doublet policy validation.status must be valid. Observed: ",
    validation_status
  )
}

if (!isTRUE(validation$validated_for_doublet_filter)) {
  stop(
    "Doublet policy has not been validated for DOUBLET_FILTER"
  )
}

if (!identical(application_status, "not_applied")) {
  stop(
    "Doublet policy application.status must be not_applied. Observed: ",
    application_status
  )
}

if (!identical(isTRUE(application$filtering_applied), FALSE)) {
  stop(
    "Doublet policy application.filtering_applied must be false"
  )
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

sample_policies <- require_mapping(
  policy$samples,
  "doublet_policy.samples"
)

if (length(sample_policies) == 0) {
  stop("doublet_policy.samples cannot be empty")
}

message("Project ID: ", policy_project_id)
message("Biological run ID: ", policy_run_id)
message("Approval status: ", approval_status)

message_step("Loading doublet-scored object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

cells_input <- ncol(seu)
features_input <- nrow(seu)
cell_ids_input <- colnames(seu)

if (cells_input == 0) {
  stop("Input Seurat object contains no cells")
}

metadata <- seu@meta.data

if (!identical(rownames(metadata), cell_ids_input)) {
  stop(
    "Metadata row names do not match Seurat cell IDs"
  )
}

required_metadata <- c(
  "doublet_class",
  "is_doublet",
  "doublet_detection_status",
  "doublet_decision_status"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required doublet metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

sample_column <- find_sample_column(metadata)
sample_values <- as.character(metadata[[sample_column]])

if (anyNA(sample_values) || any(!nzchar(sample_values))) {
  stop("Sample metadata contains missing or empty values")
}

object_samples <- sort(unique(sample_values))
policy_samples <- sort(names(sample_policies))

if (!setequal(object_samples, policy_samples)) {
  stop(
    "Policy samples do not match object samples. policy=",
    paste(policy_samples, collapse = ","),
    "; object=",
    paste(object_samples, collapse = ",")
  )
}

metadata$is_doublet <- normalise_logical(
  metadata$is_doublet
)

metadata$doublet_filter_action <- NA_character_
metadata$doublet_filter_applied <- FALSE
metadata$doublet_filter_removed <- FALSE
metadata$doublet_filter_policy_status <- approval_status

message_step("Applying per-sample doublet policy")

sample_reports <- list()
removed_reports <- list()

for (sample_name in policy_samples) {
  sample_policy <- require_mapping(
    sample_policies[[sample_name]],
    paste0("samples.", sample_name)
  )

  action <- require_scalar(
    sample_policy$application_action,
    paste0("samples.", sample_name, ".application_action")
  )

  if (!action %in% allowed_actions) {
    stop(
      "Unsupported doublet filter action for sample ",
      sample_name,
      ": ",
      action
    )
  }

  sample_index <- which(sample_values == sample_name)
  sample_cells <- cell_ids_input[sample_index]
  cells_scored <- require_int(
    sample_policy$cells_scored,
    paste0("samples.", sample_name, ".cells_scored")
  )
  approved_doublets <- require_int(
    sample_policy$proposed_doublet_count,
    paste0("samples.", sample_name, ".proposed_doublet_count")
  )

  if (length(sample_index) != cells_scored) {
    stop(
      "cells_scored mismatch for sample ",
      sample_name,
      ": policy=",
      cells_scored,
      ", object=",
      length(sample_index)
    )
  }

  detected_doublets <- sum(
    metadata$is_doublet[sample_index],
    na.rm = TRUE
  )

  if (detected_doublets != approved_doublets) {
    stop(
      "proposed_doublet_count mismatch for sample ",
      sample_name,
      ": policy=",
      approved_doublets,
      ", object=",
      detected_doublets
    )
  }

  remove_index <- rep(FALSE, length(sample_index))

  if (identical(action, "remove_proposed_doublets")) {
    remove_index <- metadata$is_doublet[sample_index]

    if (anyNA(remove_index)) {
      stop(
        "Sample ",
        sample_name,
        " has missing is_doublet values under remove_proposed_doublets"
      )
    }

    metadata$doublet_filter_applied[sample_index] <- TRUE
  }

  if (identical(action, "keep_all")) {
    metadata$doublet_filter_applied[sample_index] <- TRUE
  }

  if (identical(action, "exclude_sample_from_filtering")) {
    metadata$doublet_filter_applied[sample_index] <- FALSE
  }

  metadata$doublet_filter_action[sample_index] <- action
  metadata$doublet_filter_removed[sample_index] <- remove_index

  removed_cells <- sample_cells[remove_index]

  sample_reports[[sample_name]] <- tibble::tibble(
    sample = sample_name,
    action = action,
    cells_before = length(sample_index),
    approved_doublets = approved_doublets,
    cells_removed = length(removed_cells),
    cells_after = length(sample_index) - length(removed_cells),
    filtering_applied = any(
      metadata$doublet_filter_applied[sample_index]
    )
  )

  if (length(removed_cells) > 0) {
    removed_metadata <- metadata[removed_cells, , drop = FALSE]
    removed_reports[[sample_name]] <- tibble::tibble(
      cell_id = removed_cells,
      sample = sample_name,
      action = action,
      doublet_class = as.character(
        removed_metadata$doublet_class
      ),
      is_doublet = as.logical(
        removed_metadata$is_doublet
      ),
      doublet_detection_status = as.character(
        removed_metadata$doublet_detection_status
      ),
      doublet_decision_status = as.character(
        removed_metadata$doublet_decision_status
      )
    )
  }
}

remove_cells <- metadata$doublet_filter_removed

if (anyNA(remove_cells)) {
  stop("Internal error: missing doublet removal decisions")
}

cells_removed <- sum(remove_cells)
cells_output <- cells_input - cells_removed

expected_removed <- sum(vapply(
  policy_samples,
  function(sample_name) {
    sample_policy <- sample_policies[[sample_name]]
    action <- require_scalar(
      sample_policy$application_action,
      paste0("samples.", sample_name, ".application_action")
    )

    if (!identical(action, "remove_proposed_doublets")) {
      return(0L)
    }

    require_int(
      sample_policy$proposed_doublet_count,
      paste0("samples.", sample_name, ".proposed_doublet_count")
    )
  },
  integer(1)
))

if (cells_removed != expected_removed) {
  stop(
    "Removed cell count does not match approved removal total: ",
    cells_removed,
    " != ",
    expected_removed
  )
}

if (cells_output <= 0) {
  stop("No cells remain after applying doublet policy")
}

message("Cells before: ", cells_input)
message("Cells removed: ", cells_removed)
message("Cells after: ", cells_output)

message_step("Filtering Seurat object")

seu@meta.data <- metadata
kept_cells <- cell_ids_input[!remove_cells]

seu_filtered <- subset(
  seu,
  cells = kept_cells
)

gc()

features_output <- nrow(seu_filtered)

if (ncol(seu_filtered) != cells_output) {
  stop(
    "Filtered object cell count does not match removal decisions"
  )
}

if (features_output != features_input) {
  stop(
    "DOUBLET_FILTER unexpectedly changed the feature count"
  )
}

if (!identical(colnames(seu_filtered), kept_cells)) {
  stop(
    "Filtered object cell identity/order mismatch"
  )
}

if (cells_output != cells_input - cells_removed) {
  stop(
    "Cell count invariant failed after doublet filtering"
  )
}

message_step("Creating reports")

global_summary <- tibble::tibble(
  metric = c(
    "cells_input",
    "cells_output",
    "cells_removed",
    "features_input",
    "features_output",
    "samples",
    "approval_status",
    "validation_status",
    "application_status"
  ),
  value = c(
    cells_input,
    cells_output,
    cells_removed,
    features_input,
    features_output,
    length(policy_samples),
    approval_status,
    validation_status,
    "applied"
  )
)

by_sample <- dplyr::bind_rows(sample_reports)

removed_cells <- if (length(removed_reports) > 0) {
  dplyr::bind_rows(removed_reports)
} else {
  tibble::tibble(
    cell_id = character(),
    sample = character(),
    action = character(),
    doublet_class = character(),
    is_doublet = logical(),
    doublet_detection_status = character(),
    doublet_decision_status = character()
  )
}

message_step("Writing applied policy")

applied_at <- format(
  as.POSIXct(Sys.time(), tz = "UTC"),
  "%Y-%m-%dT%H:%M:%SZ"
)

policy$application <- list(
  status = "applied",
  filtering_applied = TRUE,
  cells_removed = cells_removed,
  applied_at = applied_at
)

policy_document$doublet_policy <- policy

message_step("Saving outputs")

saveRDS(
  seu_filtered,
  opt$output
)

readr::write_csv(
  global_summary,
  "02_doublet_filter_summary.csv"
)

readr::write_csv(
  by_sample,
  "02_doublet_filter_by_sample.csv"
)

readr::write_csv(
  removed_cells,
  "02_doublet_filter_removed_cells.csv"
)

yaml::write_yaml(
  policy_document,
  "applied_doublet_policy.yaml"
)

sink("02_doublet_filter_sessionInfo.txt")
print(sessionInfo())
sink()

message("Filtered object: ", opt$output)
message("Cells removed: ", cells_removed)
message_step("Policy-driven doublet filtering completed")
