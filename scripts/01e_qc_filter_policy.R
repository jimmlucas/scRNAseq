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
    help = "QC-metrics Seurat RDS object"
  ),
  make_option(
    "--policy",
    type = "character",
    default = NULL,
    help = "Validated and approved QC threshold policy YAML"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_qc_filtered.rds",
    help = "Output filtered Seurat RDS object"
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

normalise_optional_number <- function(
    value,
    label
) {
  if (is.null(value)) {
    return(NULL)
  }

  value <- as.numeric(value)

  if (
    length(value) != 1 ||
    is.na(value) ||
    !is.finite(value)
  ) {
    stop(label, " must be one finite number or null")
  }

  value
}

collapse_reasons <- function(reason_matrix) {
  apply(
    reason_matrix,
    1,
    function(values) {
      reasons <- names(values)[as.logical(values)]

      if (length(reasons) == 0) {
        return("")
      }

      paste(reasons, collapse = ";")
    }
  )
}

check_file(opt$input, "Input Seurat object")
check_file(opt$policy, "Approved QC policy")

message_step("Loading approved policy")

policy_document <- yaml::read_yaml(opt$policy)
policy <- require_mapping(
  policy_document$qc_threshold_policy,
  "qc_threshold_policy"
)

approval <- require_mapping(
  policy$approval,
  "qc_threshold_policy.approval"
)

validation <- require_mapping(
  policy$validation,
  "qc_threshold_policy.validation"
)

approval_status <- tolower(
  as.character(approval$status)
)

accepted_statuses <- c(
  "approved",
  "modified_and_approved"
)

if (!approval_status %in% accepted_statuses) {
  stop(
    "QC threshold policy is not approved. Status: ",
    approval_status
  )
}

if (!isTRUE(validation$validated_for_qc_filter)) {
  stop(
    "Policy has not been validated for QC_FILTER"
  )
}

if (!identical(
  tolower(as.character(validation$status)),
  "valid"
)) {
  stop(
    "Policy validation.status must be valid"
  )
}

application_semantics <- require_mapping(
  policy$application_semantics,
  "qc_threshold_policy.application_semantics"
)

if (
  !identical(
    as.character(
      application_semantics$minimum_operator
    ),
    ">="
  )
) {
  stop(
    "Unsupported minimum operator. Expected >="
  )
}

if (
  !identical(
    as.character(
      application_semantics$maximum_operator
    ),
    "<="
  )
) {
  stop(
    "Unsupported maximum operator. Expected <="
  )
}

if (
  !identical(
    toupper(
      as.character(
        application_semantics$combine_metrics_with
      )
    ),
    "AND"
  )
) {
  stop(
    "Unsupported metric combination. Expected AND"
  )
}

thresholds <- require_mapping(
  policy$thresholds,
  "qc_threshold_policy.thresholds"
)

message("Plan ID: ", policy$plan_id)
message("Project ID: ", policy$project_id)
message(
  "Biological run ID: ",
  policy$biological_run_id
)
message("Approval status: ", approval_status)
message(
  "Reviewed by: ",
  as.character(approval$reviewed_by)
)

message_step("Loading QC metrics object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

cells_before <- ncol(seu)
features_before <- nrow(seu)
cell_ids_before <- colnames(seu)

if (cells_before == 0) {
  stop("Input Seurat object contains no cells")
}

metadata <- seu@meta.data

if (!identical(
  rownames(metadata),
  cell_ids_before
)) {
  stop(
    "Metadata row names do not match Seurat cell IDs"
  )
}

enabled_thresholds <- list()

for (metric_name in names(thresholds)) {
  threshold <- require_mapping(
    thresholds[[metric_name]],
    paste0("thresholds.", metric_name)
  )

  enabled <- isTRUE(threshold$enabled)

  if (!enabled) {
    next
  }

  if (!metric_name %in% colnames(metadata)) {
    stop(
      "Enabled QC metric is missing from metadata: ",
      metric_name
    )
  }

  lower <- normalise_optional_number(
    threshold$min,
    paste0("thresholds.", metric_name, ".min")
  )

  upper <- normalise_optional_number(
    threshold$max,
    paste0("thresholds.", metric_name, ".max")
  )

  if (is.null(lower) && is.null(upper)) {
    stop(
      "Enabled threshold has no bounds: ",
      metric_name
    )
  }

  if (
    !is.null(lower) &&
    !is.null(upper) &&
    lower >= upper
  ) {
    stop(
      "Invalid bounds for ",
      metric_name,
      ": min=",
      lower,
      ", max=",
      upper
    )
  }

  metric_values <- metadata[[metric_name]]

  if (!is.numeric(metric_values)) {
    stop(
      "Enabled QC metric is not numeric: ",
      metric_name
    )
  }

  if (anyNA(metric_values)) {
    stop(
      "Enabled QC metric contains missing values: ",
      metric_name
    )
  }

  enabled_thresholds[[metric_name]] <- list(
    min = lower,
    max = upper
  )
}

if (length(enabled_thresholds) == 0) {
  stop("Approved policy has no enabled thresholds")
}

message(
  "Enabled thresholds: ",
  paste(names(enabled_thresholds), collapse = ", ")
)

message_step("Evaluating cell-level QC decisions")

reason_columns <- list()
pass_matrix <- list()

for (metric_name in names(enabled_thresholds)) {
  bounds <- enabled_thresholds[[metric_name]]
  values <- metadata[[metric_name]]

  metric_pass <- rep(TRUE, length(values))

  if (!is.null(bounds$min)) {
    reason_name <- paste0(metric_name, "_below_min")
    failed <- values < bounds$min

    reason_columns[[reason_name]] <- failed
    metric_pass <- metric_pass & !failed
  }

  if (!is.null(bounds$max)) {
    reason_name <- paste0(metric_name, "_above_max")
    failed <- values > bounds$max

    reason_columns[[reason_name]] <- failed
    metric_pass <- metric_pass & !failed
  }

  pass_matrix[[metric_name]] <- metric_pass
}

reason_df <- as.data.frame(
  reason_columns,
  check.names = FALSE
)

pass_df <- as.data.frame(
  pass_matrix,
  check.names = FALSE
)

qc_pass <- apply(
  pass_df,
  1,
  all
)

qc_fail <- !qc_pass

if (ncol(reason_df) > 0) {
  exclusion_reason <- collapse_reasons(
    reason_df
  )
} else {
  exclusion_reason <- rep(
    "",
    cells_before
  )
}

decision_metadata <- tibble::tibble(
  cell_id = cell_ids_before,
  qc_pass = qc_pass,
  qc_fail = qc_fail,
  qc_exclusion_reason = exclusion_reason
)

for (column_name in colnames(reason_df)) {
  decision_metadata[[column_name]] <-
    reason_df[[column_name]]
}

metadata$qc_pass <- qc_pass
metadata$qc_fail <- qc_fail
metadata$qc_exclusion_reason <-
  exclusion_reason

seu@meta.data <- metadata

cells_removed <- sum(qc_fail)
cells_after <- sum(qc_pass)

message("Cells before: ", cells_before)
message("Cells passing: ", cells_after)
message("Cells removed: ", cells_removed)

if (cells_after == 0) {
  stop(
    "No cells remain after applying approved QC policy"
  )
}

message_step("Filtering Seurat object")

filtered_cells <- cell_ids_before[qc_pass]

seu_filtered <- subset(
  seu,
  cells = filtered_cells
)

gc()

if (ncol(seu_filtered) != cells_after) {
  stop(
    "Filtered object cell count does not match QC decisions"
  )
}

if (nrow(seu_filtered) != features_before) {
  stop(
    "QC_FILTER unexpectedly changed the feature count"
  )
}

if (!identical(
  colnames(seu_filtered),
  filtered_cells
)) {
  stop(
    "Filtered object cell identity/order mismatch"
  )
}

message_step("Creating reports")

global_summary <- tibble::tibble(
  metric = c(
    "cells_before_filtering",
    "cells_after_filtering",
    "cells_removed",
    "fraction_removed",
    "features_before_filtering",
    "features_after_filtering",
    "enabled_threshold_count"
  ),
  value = c(
    cells_before,
    cells_after,
    cells_removed,
    cells_removed / cells_before,
    features_before,
    nrow(seu_filtered),
    length(enabled_thresholds)
  )
)

threshold_rows <- lapply(
  names(enabled_thresholds),
  function(metric_name) {
    bounds <- enabled_thresholds[[metric_name]]
    values <- metadata[[metric_name]]

    below_count <- if (is.null(bounds$min)) {
      0L
    } else {
      sum(values < bounds$min)
    }

    above_count <- if (is.null(bounds$max)) {
      0L
    } else {
      sum(values > bounds$max)
    }

    tibble::tibble(
      metric = metric_name,
      min = if (is.null(bounds$min)) {
        NA_real_
      } else {
        bounds$min
      },
      max = if (is.null(bounds$max)) {
        NA_real_
      } else {
        bounds$max
      },
      cells_below_min = below_count,
      cells_above_max = above_count,
      cells_flagged_for_metric = sum(
        !pass_df[[metric_name]]
      )
    )
  }
)

threshold_summary <- dplyr::bind_rows(
  threshold_rows
)

sample_column <- if (
  "sample" %in% colnames(metadata)
) {
  "sample"
} else {
  NULL
}

group_column <- if (
  "Group" %in% colnames(metadata)
) {
  "Group"
} else {
  NULL
}

sample_report <- metadata

if (!"cell_id" %in% colnames(sample_report)) {
  sample_report$cell_id <- rownames(sample_report)
}

if (!is.null(sample_column)) {
  grouping_columns <- sample_column

  if (!is.null(group_column)) {
    grouping_columns <- c(
      grouping_columns,
      group_column
    )
  }

  sample_report <- sample_report %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(grouping_columns)
      )
    ) %>%
    dplyr::summarise(
      cells_before = dplyr::n(),
      cells_after = sum(qc_pass),
      cells_removed = sum(qc_fail),
      fraction_removed = mean(qc_fail),
      .groups = "drop"
    )
} else {
  sample_report <- tibble::tibble(
    sample = "all_cells",
    cells_before = cells_before,
    cells_after = cells_after,
    cells_removed = cells_removed,
    fraction_removed = cells_removed / cells_before
  )
}

message_step("Saving outputs")

saveRDS(
  seu_filtered,
  opt$output
)

readr::write_csv(
  global_summary,
  "01e_qc_filter_summary.csv"
)

readr::write_csv(
  threshold_summary,
  "01e_qc_filter_by_threshold.csv"
)

readr::write_csv(
  sample_report,
  "01e_qc_filter_by_sample.csv"
)

readr::write_csv(
  decision_metadata,
  "01e_qc_filter_cell_decisions.csv"
)

sink("01e_sessionInfo.txt")
print(sessionInfo())
sink()

message("Filtered object: ", opt$output)
message(
  "Cell decision report: ",
  "01e_qc_filter_cell_decisions.csv"
)
message_step(
  "Policy-driven QC filtering completed"
)
