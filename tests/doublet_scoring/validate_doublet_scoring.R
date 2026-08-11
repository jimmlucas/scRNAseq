#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  "tests/doublet_scoring/fixtures/seu_qc_filtered.rds"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  "results_test_repro/01g_doublet_detection_scoring"
}

output_rds <- file.path(
  output_dir,
  "seu_doublet_scored.rds"
)

summary_csv <- file.path(
  output_dir,
  "01g_doublet_detection_summary.csv"
)

by_sample_csv <- file.path(
  output_dir,
  "01g_doublet_detection_by_sample.csv"
)

cell_scores_csv <- file.path(
  output_dir,
  "01g_doublet_detection_cell_scores.csv"
)

performance_csv <- file.path(
  output_dir,
  "01g_doublet_detection_performance.csv"
)

pk_sweep_csv <- file.path(
  output_dir,
  "01g_doublet_detection_pk_sweep.csv"
)

required_files <- c(
  input_path,
  output_rds,
  summary_csv,
  by_sample_csv,
  cell_scores_csv,
  performance_csv,
  pk_sweep_csv,
  file.path(
    output_dir,
    "01g_sessionInfo.txt"
  )
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

before <- readRDS(input_path)
after <- readRDS(output_rds)

if (!inherits(after, "Seurat")) {
  stop("Output object is not a Seurat object.")
}

# Invariant preservation contract.
stopifnot(
  ncol(before) == 392,
  nrow(before) == 27522,
  ncol(after) == ncol(before),
  nrow(after) == nrow(before),
  identical(colnames(before), colnames(after)),
  identical(rownames(before), rownames(after))
)

required_metadata <- c(
  "doublet_score",
  "doublet_class",
  "is_doublet",
  "doublet_pK",
  "doublet_nExp_raw",
  "doublet_nExp_adjusted",
  "doublet_detection_status",
  "doublet_decision_status",
  "doublet_reviewed_class"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(after@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

detection_status <- table(
  after$doublet_detection_status,
  useNA = "ifany"
)

decision_status <- table(
  after$doublet_decision_status,
  useNA = "ifany"
)

doublet_classes <- table(
  after$doublet_class,
  useNA = "ifany"
)

stopifnot(
  identical(
    unname(detection_status[["completed"]]),
    392L
  ),
  identical(
    unname(decision_status[["unreviewed"]]),
    392L
  ),
  identical(
    unname(doublet_classes[["Doublet"]]),
    29L
  ),
  identical(
    unname(doublet_classes[["Singlet"]]),
    363L
  )
)

summary_table <- read.csv(
  summary_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_values <- setNames(
  summary_table$value,
  summary_table$metric
)

stopifnot(
  as.integer(summary_values[["cells_input"]]) == 392L,
  as.integer(summary_values[["cells_output"]]) == 392L,
  as.integer(summary_values[["features_input"]]) == 27522L,
  as.integer(summary_values[["features_output"]]) == 27522L,
  as.integer(summary_values[["cells_removed"]]) == 0L,
  as.integer(summary_values[["samples"]]) == 2L,
  as.integer(summary_values[["proposed_doublets"]]) == 29L
)

by_sample <- read.csv(
  by_sample_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(by_sample) == 2L,
  setequal(by_sample$sample, c("wta-1", "wta-5")),
  all(
    by_sample$doublet_detection_status == "completed"
  ),
  all(
    by_sample$doublet_decision_status == "unreviewed"
  ),
  sum(by_sample$cells_scored) == 392L,
  sum(by_sample$proposed_doublets) == 29L,
  all(by_sample$pk_mode == "auto"),
  all(by_sample$pk_source == "adaptive_stable"),
  all(by_sample$pk_stable),
  all(by_sample$sweep_cells == by_sample$cells_scored),
  all(by_sample$sweep_attempts == 2L)
)

pk_sweep <- read.csv(
  pk_sweep_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  all(pk_sweep$pk_mode == "auto"),
  all(pk_sweep$sampling_method == "all_cells"),
  all(pk_sweep$sweep_seed %in% c(1234L, 4321L)),
  sum(pk_sweep$selected_final) == 2L
)

cell_scores <- read.csv(
  cell_scores_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(cell_scores) == 392L,
  sum(cell_scores$doublet_class == "Doublet") == 29L,
  sum(cell_scores$doublet_class == "Singlet") == 363L,
  all(
    cell_scores$doublet_detection_status == "completed"
  ),
  all(
    cell_scores$doublet_decision_status == "unreviewed"
  )
)

performance <- read.csv(
  performance_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_performance_columns <- c(
  "sample",
  "cells",
  "nfeatures",
  "pcs",
  "sweep_cores",
  "pk_mode",
  "pk_source",
  "pk_stable",
  "pk_delta",
  "selected_bcmetric",
  "sweep_cells",
  "sweep_fraction",
  "sampling_method",
  "sweep_attempts",
  "normalize_seconds",
  "find_variable_features_seconds",
  "scale_data_seconds",
  "run_pca_seconds",
  "param_sweep_seconds",
  "doublet_finder_seconds",
  "total_seconds"
)

missing_performance_columns <- setdiff(
  required_performance_columns,
  colnames(performance)
)

if (length(missing_performance_columns) > 0) {
  stop(
    "Missing required performance columns: ",
    paste(missing_performance_columns, collapse = ", ")
  )
}

stopifnot(
  nrow(performance) == 2L,
  setequal(performance$sample, c("wta-1", "wta-5")),
  sum(performance$cells) == 392L,
  all(performance$nfeatures <= 1000L),
  all(performance$pcs == 10L),
  all(performance$sweep_cores >= 1L),
  all(performance$pk_mode == "auto"),
  all(performance$pk_stable),
  all(performance$sweep_cells == performance$cells),
  all(performance$sweep_attempts == 2L),
  all(performance$total_seconds >= 0),
  all(performance$param_sweep_seconds >= 0)
)

cat("\nDoublet scoring regression test: VALID\n")
cat("Cells preserved: 392\n")
cat("Features preserved: 27522\n")
cat("Doublets proposed: 29\n")
cat("Singlets proposed: 363\n")
cat("Cells removed: 0\n")
cat("Performance metrics: recorded\n")
