#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option(
    "--before-object",
    dest = "before_object",
    type = "character",
    default = NULL,
    help = "Seurat object containing QC metrics before filtering"
  ),
  make_option(
    "--after-object",
    dest = "after_object",
    type = "character",
    default = NULL,
    help = "Filtered Seurat object"
  ),
  make_option(
    "--cell-decisions",
    dest = "cell_decisions",
    type = "character",
    default = NULL,
    help = "Cell-level QC filtering decisions CSV"
  ),
  make_option(
    "--policy",
    type = "character",
    default = NULL,
    help = "Validated approved QC threshold policy YAML"
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

safe_save_plot <- function(
    plot,
    filename,
    width = 8,
    height = 5,
    dpi = 150
) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  message("Plot saved: ", filename)
}

normalise_metadata <- function(
    seu,
    stage_label
) {
  metadata <- seu@meta.data

  if (!"cell_id" %in% colnames(metadata)) {
    metadata$cell_id <- rownames(metadata)
  }

  if (!"sample" %in% colnames(metadata)) {
    if ("ShortName" %in% colnames(metadata)) {
      metadata$sample <- metadata$ShortName
    } else if ("sample_id" %in% colnames(metadata)) {
      metadata$sample <- metadata$sample_id
    } else if ("orig.ident" %in% colnames(metadata)) {
      metadata$sample <- metadata$orig.ident
    } else {
      metadata$sample <- "sample"
    }
  }

  if (!"Group" %in% colnames(metadata)) {
    metadata$Group <- NA_character_
  }

  metadata$qc_stage <- stage_label
  metadata
}

summarise_metric <- function(
    data,
    metric_name
) {
  values <- data[[metric_name]]

  tibble::tibble(
    qc_stage = unique(data$qc_stage),
    metric = metric_name,
    cells = length(values),
    minimum = min(values, na.rm = TRUE),
    q25 = as.numeric(
      stats::quantile(
        values,
        probs = 0.25,
        na.rm = TRUE
      )
    ),
    median = median(values, na.rm = TRUE),
    mean = mean(values, na.rm = TRUE),
    q75 = as.numeric(
      stats::quantile(
        values,
        probs = 0.75,
        na.rm = TRUE
      )
    ),
    maximum = max(values, na.rm = TRUE)
  )
}

check_file(
  opt$before_object,
  "Before-filtering Seurat object"
)
check_file(
  opt$after_object,
  "After-filtering Seurat object"
)
check_file(
  opt$cell_decisions,
  "Cell decisions report"
)
check_file(
  opt$policy,
  "Approved QC policy"
)

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

if (!approval_status %in% c(
  "approved",
  "modified_and_approved"
)) {
  stop(
    "QC visualization requires an approved policy. Status: ",
    approval_status
  )
}

if (!isTRUE(validation$validated_for_qc_filter)) {
  stop(
    "Policy has not been validated for QC_FILTER"
  )
}

message_step("Loading traced QC artifacts")

before_seu <- readRDS(opt$before_object)
after_seu <- readRDS(opt$after_object)

if (!inherits(before_seu, "Seurat")) {
  stop("Before-filtering object is not a Seurat object")
}

if (!inherits(after_seu, "Seurat")) {
  stop("After-filtering object is not a Seurat object")
}

before_cells <- colnames(before_seu)
after_cells <- colnames(after_seu)

if (!all(after_cells %in% before_cells)) {
  stop(
    "After-filtering cells are not a subset of before-filtering cells"
  )
}

if (nrow(before_seu) != nrow(after_seu)) {
  stop(
    "Feature count differs between before and after objects"
  )
}

decisions <- readr::read_csv(
  opt$cell_decisions,
  show_col_types = FALSE
)

required_decision_columns <- c(
  "cell_id",
  "qc_pass",
  "qc_fail",
  "qc_exclusion_reason"
)

missing_decision_columns <- setdiff(
  required_decision_columns,
  colnames(decisions)
)

if (length(missing_decision_columns) > 0) {
  stop(
    "Missing decision columns: ",
    paste(
      missing_decision_columns,
      collapse = ", "
    )
  )
}

if (!setequal(decisions$cell_id, before_cells)) {
  stop(
    "Cell decision report does not match the before-filtering object"
  )
}

expected_after_cells <- decisions %>%
  dplyr::filter(qc_pass) %>%
  dplyr::pull(cell_id)

if (!identical(after_cells, expected_after_cells)) {
  stop(
    "After-filtering object does not match the approved cell decisions"
  )
}

before_meta <- normalise_metadata(
  before_seu,
  "before"
)

after_meta <- normalise_metadata(
  after_seu,
  "after"
)

required_metrics <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mito",
  "percent.ribo",
  "percent.globin",
  "log10GenesPerUMI"
)

for (metric_name in required_metrics) {
  if (!metric_name %in% colnames(before_meta)) {
    stop(
      "Before-filtering metadata is missing metric: ",
      metric_name
    )
  }

  if (!metric_name %in% colnames(after_meta)) {
    stop(
      "After-filtering metadata is missing metric: ",
      metric_name
    )
  }
}

message_step("Creating traceable summaries")

combined_metadata <- dplyr::bind_rows(
  before_meta,
  after_meta
)

metric_summary <- dplyr::bind_rows(
  lapply(
    required_metrics,
    function(metric_name) {
      dplyr::bind_rows(
        summarise_metric(
          before_meta,
          metric_name
        ),
        summarise_metric(
          after_meta,
          metric_name
        )
      )
    }
  )
)

cell_summary <- tibble::tibble(
  metric = c(
    "cells_before",
    "cells_after",
    "cells_removed",
    "fraction_removed",
    "features_before",
    "features_after",
    "approved_threshold_count"
  ),
  value = c(
    ncol(before_seu),
    ncol(after_seu),
    ncol(before_seu) - ncol(after_seu),
    (
      ncol(before_seu) - ncol(after_seu)
    ) / ncol(before_seu),
    nrow(before_seu),
    nrow(after_seu),
    as.numeric(
      validation$enabled_threshold_count
    )
  )
)

sample_summary <- combined_metadata %>%
  dplyr::group_by(
    qc_stage,
    sample,
    Group
  ) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nCount_RNA = median(
      nCount_RNA,
      na.rm = TRUE
    ),
    median_nFeature_RNA = median(
      nFeature_RNA,
      na.rm = TRUE
    ),
    median_percent_mito = median(
      percent.mito,
      na.rm = TRUE
    ),
    median_percent_ribo = median(
      percent.ribo,
      na.rm = TRUE
    ),
    median_percent_globin = median(
      percent.globin,
      na.rm = TRUE
    ),
    median_log10GenesPerUMI = median(
      log10GenesPerUMI,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

exclusion_summary <- decisions %>%
  dplyr::filter(qc_fail) %>%
  dplyr::count(
    qc_exclusion_reason,
    name = "cells"
  ) %>%
  dplyr::arrange(
    dplyr::desc(cells)
  )

if (nrow(exclusion_summary) == 0) {
  exclusion_summary <- tibble::tibble(
    qc_exclusion_reason = character(),
    cells = integer()
  )
}

readr::write_csv(
  cell_summary,
  "01f_qc_visualization_summary.csv"
)

readr::write_csv(
  metric_summary,
  "01f_qc_metric_distributions.csv"
)

readr::write_csv(
  sample_summary,
  "01f_qc_by_sample_before_after.csv"
)

readr::write_csv(
  exclusion_summary,
  "01f_qc_exclusion_reasons.csv"
)

readr::write_csv(
  combined_metadata,
  "01f_qc_cell_metadata_before_after.csv"
)

message_step("Generating before/after QC plots")

metric_labels <- c(
  nCount_RNA = "UMI count",
  nFeature_RNA = "Detected genes",
  percent.mito = "Mitochondrial percentage",
  percent.ribo = "Ribosomal percentage",
  percent.globin = "Globin percentage",
  log10GenesPerUMI = "Library complexity"
)

for (metric_name in required_metrics) {
  plot_data <- combined_metadata %>%
    dplyr::select(
      qc_stage,
      sample,
      Group,
      dplyr::all_of(metric_name)
    )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = sample,
      y = .data[[metric_name]],
      fill = qc_stage
    )
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      position = ggplot2::position_dodge(
        width = 0.8
      )
    ) +
    ggplot2::geom_boxplot(
      width = 0.12,
      outlier.size = 0.2,
      position = ggplot2::position_dodge(
        width = 0.8
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    ) +
    ggplot2::labs(
      title = paste0(
        metric_labels[[metric_name]],
        " before and after QC filtering"
      ),
      x = "Sample",
      y = metric_labels[[metric_name]],
      fill = "QC stage"
    )

  safe_save_plot(
    plot_object,
    paste0(
      "01f_qc_",
      gsub("\\.", "_", metric_name),
      "_before_after.png"
    )
  )
}

message_step("Generating sample retention plot")

sample_retention <- decisions %>%
  dplyr::left_join(
    before_meta %>%
      dplyr::select(
        cell_id,
        sample,
        Group
      ),
    by = "cell_id"
  ) %>%
  dplyr::group_by(
    sample,
    Group
  ) %>%
  dplyr::summarise(
    cells_before = dplyr::n(),
    cells_after = sum(qc_pass),
    cells_removed = sum(qc_fail),
    fraction_retained = mean(qc_pass),
    .groups = "drop"
  )

readr::write_csv(
  sample_retention,
  "01f_qc_sample_retention.csv"
)

retention_long <- sample_retention %>%
  dplyr::select(
    sample,
    Group,
    cells_before,
    cells_after
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      cells_before,
      cells_after
    ),
    names_to = "qc_stage",
    values_to = "cells"
  ) %>%
  dplyr::mutate(
    qc_stage = dplyr::recode(
      qc_stage,
      cells_before = "before",
      cells_after = "after"
    )
  )

p_retention <- ggplot2::ggplot(
  retention_long,
  ggplot2::aes(
    x = sample,
    y = cells,
    fill = qc_stage
  )
) +
  ggplot2::geom_col(
    position = "dodge"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Cell retention after approved QC filtering",
    x = "Sample",
    y = "Cell count",
    fill = "QC stage"
  )

safe_save_plot(
  p_retention,
  "01f_qc_cell_retention_by_sample.png"
)

if (nrow(exclusion_summary) > 0) {
  p_reasons <- ggplot2::ggplot(
    exclusion_summary,
    ggplot2::aes(
      x = reorder(
        qc_exclusion_reason,
        cells
      ),
      y = cells
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "QC exclusion reasons",
      x = "Exclusion reason",
      y = "Cells"
    )

  safe_save_plot(
    p_reasons,
    "01f_qc_exclusion_reasons.png"
  )
}

message_step("Generating post-filter scatter plot")

p_scatter <- ggplot2::ggplot(
  after_meta,
  ggplot2::aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    colour = percent.mito
  )
) +
  ggplot2::geom_point(
    size = 0.5,
    alpha = 0.6
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Cell complexity after approved QC filtering",
    x = "nCount_RNA",
    y = "nFeature_RNA",
    colour = "percent.mito"
  )

safe_save_plot(
  p_scatter,
  "01f_qc_scatter_after_filtering.png",
  width = 7,
  height = 5
)

sink("01f_sessionInfo.txt")
print(sessionInfo())
sink()

message("Cells before: ", ncol(before_seu))
message("Cells after: ", ncol(after_seu))
message(
  "Approval status: ",
  approval_status
)
message_step(
  "Traceable QC visualization completed"
)
