#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  "results_test_repro/01g_doublet_detection_scoring/seu_doublet_scored.rds"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  "results_test_repro/02_doublet_filter"
}

output_rds <- file.path(
  output_dir,
  "seu_doublet_filtered.rds"
)

summary_csv <- file.path(
  output_dir,
  "02_doublet_filter_summary.csv"
)

by_sample_csv <- file.path(
  output_dir,
  "02_doublet_filter_by_sample.csv"
)

removed_csv <- file.path(
  output_dir,
  "02_doublet_filter_removed_cells.csv"
)

applied_policy <- file.path(
  output_dir,
  "applied_doublet_policy.yaml"
)

required_files <- c(
  input_path,
  output_rds,
  summary_csv,
  by_sample_csv,
  removed_csv,
  applied_policy,
  file.path(
    output_dir,
    "02_doublet_filter_sessionInfo.txt"
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

stopifnot(
  ncol(before) == 392L,
  nrow(before) == 27522L,
  ncol(after) == 363L,
  nrow(after) == nrow(before),
  identical(rownames(before), rownames(after))
)

required_metadata <- c(
  "doublet_class",
  "is_doublet",
  "doublet_detection_status",
  "doublet_decision_status",
  "doublet_filter_action",
  "doublet_filter_applied",
  "doublet_filter_removed",
  "doublet_filter_policy_status"
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

stopifnot(
  all(after$doublet_filter_removed == FALSE),
  all(after$doublet_filter_applied == TRUE),
  all(after$doublet_filter_action == "remove_proposed_doublets"),
  all(after$doublet_filter_policy_status == "approved"),
  all(after$is_doublet == FALSE)
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
  as.integer(summary_values[["cells_output"]]) == 363L,
  as.integer(summary_values[["cells_removed"]]) == 29L,
  as.integer(summary_values[["features_input"]]) == 27522L,
  as.integer(summary_values[["features_output"]]) == 27522L,
  summary_values[["approval_status"]] == "approved",
  summary_values[["validation_status"]] == "valid",
  summary_values[["application_status"]] == "applied"
)

by_sample <- read.csv(
  by_sample_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(by_sample) == 2L,
  setequal(by_sample$sample, c("wta-1", "wta-5")),
  all(by_sample$action == "remove_proposed_doublets"),
  sum(by_sample$cells_before) == 392L,
  sum(by_sample$approved_doublets) == 29L,
  sum(by_sample$cells_removed) == 29L,
  sum(by_sample$cells_after) == 363L
)

removed <- read.csv(
  removed_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(removed) == 29L,
  all(removed$is_doublet),
  all(removed$doublet_class == "Doublet"),
  setequal(removed$sample, c("wta-1", "wta-5"))
)

applied <- yaml::read_yaml(applied_policy)
application <- applied$doublet_policy$application

stopifnot(
  applied$doublet_policy$approval$status == "approved",
  applied$doublet_policy$validation$status == "valid",
  application$status == "applied",
  isTRUE(application$filtering_applied),
  as.integer(application$cells_removed) == 29L
)

cat("\nDoublet filter regression test: VALID\n")
cat("Cells input: 392\n")
cat("Cells removed: 29\n")
cat("Cells output: 363\n")
cat("Features preserved: 27522\n")
