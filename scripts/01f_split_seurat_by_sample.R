#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option(
    "--input",
    type = "character",
    default = NULL,
    help = "QC-filtered Seurat RDS object"
  ),
  make_option(
    "--output-dir",
    dest = "output_dir",
    type = "character",
    default = "sample_objects",
    help = "Directory where per-sample Seurat objects will be written"
  ),
  make_option(
    "--manifest",
    type = "character",
    default = "01f_split_by_sample_manifest.tsv",
    help = "Output TSV manifest with one row per sample object"
  ),
  make_option(
    "--summary",
    type = "character",
    default = "01f_split_by_sample_summary.csv",
    help = "Output CSV summary"
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

safe_id <- function(value) {
  value <- as.character(value)
  value <- trimws(value)
  value <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  value <- gsub("^_+|_+$", "", value)

  if (!nzchar(value)) {
    return("sample")
  }

  value
}

detect_sample_column <- function(metadata) {
  candidates <- c(
    "sample",
    "ShortName",
    "sample_id",
    "orig.ident"
  )

  available <- candidates[candidates %in% colnames(metadata)]

  if (length(available) == 0) {
    stop(
      "Could not find a sample column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }

  available[[1]]
}

write_session_info <- function(path) {
  sink(path)
  on.exit(sink(), add = TRUE)
  print(sessionInfo())
}

check_file(opt$input, "Input Seurat object")

message_step("Loading QC-filtered Seurat object")

seu <- readRDS(opt$input)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object")
}

metadata <- seu[[]]
sample_column <- detect_sample_column(metadata)
sample_values <- as.character(metadata[[sample_column]])

if (any(is.na(sample_values)) || any(!nzchar(trimws(sample_values)))) {
  stop("Sample column contains missing or empty values: ", sample_column)
}

samples <- sort(unique(sample_values))
sample_safe_ids <- vapply(samples, safe_id, character(1))

if (anyDuplicated(sample_safe_ids)) {
  duplicated_ids <- unique(sample_safe_ids[duplicated(sample_safe_ids)])
  stop(
    "Sample ids collide after filesystem sanitization: ",
    paste(duplicated_ids, collapse = ", ")
  )
}

dir.create(
  opt$output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cells_input <- ncol(seu)
features_input <- nrow(seu)
manifest_rows <- vector("list", length(samples))

message_step("Writing per-sample Seurat objects")

for (idx in seq_along(samples)) {
  sample_id <- samples[[idx]]
  sample_safe <- sample_safe_ids[[idx]]
  cells <- rownames(metadata)[sample_values == sample_id]

  if (length(cells) == 0) {
    stop("Sample has no cells after split: ", sample_id)
  }

  sample_dir <- file.path(opt$output_dir, sample_safe)
  dir.create(
    sample_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  object_path <- file.path(
    sample_dir,
    paste0("seu_qc_filtered.", sample_safe, ".rds")
  )

  message(
    "[", sample_id, "] writing ",
    length(cells), " cells to ", object_path
  )

  sample_object <- subset(
    seu,
    cells = cells
  )

  if (nrow(sample_object) != features_input) {
    stop(
      "Feature count changed after splitting sample ",
      sample_id,
      ": expected ",
      features_input,
      ", observed ",
      nrow(sample_object)
    )
  }

  saveRDS(
    sample_object,
    object_path
  )

  manifest_rows[[idx]] <- data.frame(
    sample = sample_id,
    sample_safe_id = sample_safe,
    sample_column = sample_column,
    cells = ncol(sample_object),
    features = nrow(sample_object),
    object_path = object_path,
    stringsAsFactors = FALSE
  )

  rm(sample_object)
  gc(verbose = FALSE)
}

manifest <- do.call(
  rbind,
  manifest_rows
)

cells_output <- sum(manifest$cells)

if (cells_output != cells_input) {
  stop(
    "Split cell count mismatch: input=",
    cells_input,
    ", output=",
    cells_output
  )
}

if (any(manifest$features != features_input)) {
  stop("At least one sample object changed the feature count")
}

summary_df <- data.frame(
  metric = c(
    "cells_input",
    "features_input",
    "samples",
    "cells_output",
    "features_output",
    "sample_column"
  ),
  value = c(
    as.character(cells_input),
    as.character(features_input),
    as.character(length(samples)),
    as.character(cells_output),
    as.character(features_input),
    sample_column
  ),
  stringsAsFactors = FALSE
)

message_step("Writing manifest and summary")

write.table(
  manifest,
  opt$manifest,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  summary_df,
  opt$summary,
  row.names = FALSE,
  quote = TRUE
)

write_session_info("01f_sessionInfo.txt")

message_step("Split complete")
message(
  "Samples: ",
  length(samples),
  "; cells: ",
  cells_input,
  "; features: ",
  features_input
)
