#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  "results_test_repro/03_post_doublet_qc"
}

summary_csv <- file.path(
  output_dir,
  "03_post_doublet_qc_summary.csv"
)
by_sample_csv <- file.path(
  output_dir,
  "03_post_doublet_qc_by_sample.csv"
)
warnings_tsv <- file.path(
  output_dir,
  "03_post_doublet_qc_warnings.tsv"
)
readiness_yaml <- file.path(
  output_dir,
  "03_analysis_readiness.yaml"
)

required_files <- c(
  summary_csv,
  by_sample_csv,
  warnings_tsv,
  readiness_yaml,
  file.path(output_dir, "03_sessionInfo.txt")
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
  as.integer(summary_values[["cells"]]) == 363L,
  as.integer(summary_values[["features"]]) == 27522L,
  as.integer(summary_values[["samples"]]) == 2L,
  as.integer(summary_values[["doublets_removed"]]) == 29L,
  as.integer(summary_values[["warning_count"]]) == 0L
)

required_metrics <- c(
  "nFeature_RNA_median",
  "nCount_RNA_median",
  "percent_mito_median"
)

missing_metrics <- setdiff(
  required_metrics,
  summary_table$metric
)

if (length(missing_metrics) > 0) {
  stop(
    "Missing expected post-doublet QC metrics: ",
    paste(missing_metrics, collapse = ", ")
  )
}

by_sample <- read.csv(
  by_sample_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(by_sample) == 2L,
  setequal(by_sample$sample, c("wta-1", "wta-5")),
  all(by_sample$action == "remove_proposed_doublets"),
  sum(by_sample$cells_original) == 392L,
  sum(by_sample$doublets_removed) == 29L,
  sum(by_sample$cells_final) == 363L,
  all(by_sample$cells_final >= 50L)
)

warnings <- read.delim(
  warnings_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(warnings) == 0L)

readiness <- yaml::read_yaml(readiness_yaml)
contract <- readiness$analysis_readiness

stopifnot(
  contract$status == "ready",
  isTRUE(contract$qc$completed),
  isTRUE(contract$doublets$scoring_completed),
  isTRUE(contract$doublets$human_approval_validated),
  isTRUE(contract$doublets$filtering_applied),
  as.integer(contract$dataset$cells) == 363L,
  as.integer(contract$dataset$features) == 27522L,
  as.integer(contract$dataset$samples) == 2L,
  contract$next_phase$recommended ==
    "normalization_policy_resolution"
)

cat("\nPost-doublet QC readiness test: VALID\n")
cat("Readiness: ready\n")
cat("Cells: 363\n")
cat("Features: 27522\n")
cat("Samples: 2\n")
