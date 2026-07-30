#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(
    "--scored-pattern",
    type = "character",
    default = "seu_doublet_scored.*.rds",
    help = "Glob pattern for per-sample doublet-scored RDS objects"
  ),
  make_option(
    "--by-sample-pattern",
    type = "character",
    default = "01g_doublet_detection_by_sample.*.csv",
    help = "Glob pattern for per-sample by-sample summaries"
  ),
  make_option(
    "--cell-scores-pattern",
    type = "character",
    default = "01g_doublet_detection_cell_scores.*.csv",
    help = "Glob pattern for per-sample cell score tables"
  ),
  make_option(
    "--pk-sweep-pattern",
    type = "character",
    default = "01g_doublet_detection_pk_sweep.*.csv",
    help = "Glob pattern for per-sample pK sweep tables"
  ),
  make_option(
    "--performance-pattern",
    type = "character",
    default = "01g_doublet_detection_performance.*.csv",
    help = "Glob pattern for per-sample performance tables"
  ),
  make_option(
    "--output",
    type = "character",
    default = "seu_doublet_scored.rds",
    help = "Merged output Seurat object"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

message_step <- function(label) {
  message("\n========== ", label, " ==========")
}

list_required <- function(pattern, label) {
  paths <- Sys.glob(pattern)
  paths <- sort(paths)

  if (length(paths) == 0) {
    stop("No files found for ", label, " with pattern: ", pattern)
  }

  paths
}

read_csv_many <- function(paths) {
  dplyr::bind_rows(
    lapply(
      paths,
      function(path) {
        readr::read_csv(
          path,
          show_col_types = FALSE
        )
      }
    )
  )
}

message_step("Loading per-sample doublet-scored objects")

scored_paths <- list_required(
  opt$scored_pattern,
  "scored objects"
)

objects <- lapply(scored_paths, readRDS)

if (!all(vapply(objects, inherits, logical(1), what = "Seurat"))) {
  stop("At least one scored object is not a Seurat object")
}

feature_counts <- vapply(objects, nrow, integer(1))

if (length(unique(feature_counts)) != 1) {
  stop(
    "Per-sample scored objects do not have the same feature count: ",
    paste(feature_counts, collapse = ", ")
  )
}

message_step("Merging Seurat objects")

merged <- objects[[1]]

if (length(objects) > 1) {
  merged <- merge(
    x = objects[[1]],
    y = objects[-1],
    merge.data = TRUE
  )
}

rm(objects)
gc(verbose = FALSE)

required_metadata <- c(
  "doublet_score",
  "doublet_class",
  "is_doublet",
  "doublet_detection_status",
  "doublet_decision_status",
  "doublet_reviewed_class",
  "doublet_pK",
  "doublet_nExp_raw",
  "doublet_nExp_adjusted"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(merged[[]])
)

if (length(missing_metadata) > 0) {
  stop(
    "Merged object is missing doublet metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

message_step("Merging tabular artifacts")

by_sample <- read_csv_many(
  list_required(
    opt$by_sample_pattern,
    "by-sample summaries"
  )
)

cell_scores <- read_csv_many(
  list_required(
    opt$cell_scores_pattern,
    "cell score tables"
  )
)

pk_sweep <- read_csv_many(
  list_required(
    opt$pk_sweep_pattern,
    "pK sweep tables"
  )
)

performance <- read_csv_many(
  list_required(
    opt$performance_pattern,
    "performance tables"
  )
)

if (nrow(cell_scores) != ncol(merged)) {
  stop(
    "Cell score row count does not match merged object cells: ",
    nrow(cell_scores),
    " vs ",
    ncol(merged)
  )
}

proposed_doublets <- sum(
  as.integer(by_sample$proposed_doublets),
  na.rm = TRUE
)

cells_input <- ncol(merged)

overall_summary <- data.frame(
  metric = c(
    "cells_input",
    "cells_output",
    "features_input",
    "features_output",
    "cells_removed",
    "samples",
    "proposed_doublets",
    "proposed_doublet_fraction",
    "approval_status"
  ),
  value = c(
    cells_input,
    ncol(merged),
    nrow(merged),
    nrow(merged),
    0,
    nrow(by_sample),
    proposed_doublets,
    proposed_doublets / cells_input,
    "proposed"
  )
)

message_step("Saving merged doublet scoring outputs")

saveRDS(
  merged,
  opt$output
)

readr::write_csv(
  overall_summary,
  "01g_doublet_detection_summary.csv"
)

readr::write_csv(
  by_sample,
  "01g_doublet_detection_by_sample.csv"
)

readr::write_csv(
  cell_scores,
  "01g_doublet_detection_cell_scores.csv"
)

readr::write_csv(
  pk_sweep,
  "01g_doublet_detection_pk_sweep.csv"
)

readr::write_csv(
  performance,
  "01g_doublet_detection_performance.csv"
)

sink("01g_sessionInfo.txt")
print(sessionInfo())
sink()

message("Merged doublet-scored object: ", opt$output)
message_step("01g_merge_doublet_detection_scoring finished successfully")
