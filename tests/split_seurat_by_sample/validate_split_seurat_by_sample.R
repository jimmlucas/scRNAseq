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
    help = "Original QC-filtered Seurat RDS object"
  ),
  make_option(
    "--manifest",
    type = "character",
    default = NULL,
    help = "Split manifest TSV"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("Input object does not exist")
}

if (is.null(opt$manifest) || !file.exists(opt$manifest)) {
  stop("Manifest does not exist")
}

seu <- readRDS(opt$input)
manifest <- read.delim(
  opt$manifest,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample",
  "sample_safe_id",
  "sample_column",
  "cells",
  "features",
  "object_path"
)

missing_columns <- setdiff(
  required_columns,
  colnames(manifest)
)

if (length(missing_columns) > 0) {
  stop(
    "Manifest is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (sum(as.integer(manifest$cells)) != ncol(seu)) {
  stop("Manifest cell counts do not match input object")
}

if (any(as.integer(manifest$features) != nrow(seu))) {
  stop("Manifest feature counts do not match input object")
}

sample_column <- manifest$sample_column[[1]]

if (!sample_column %in% colnames(seu[[]])) {
  stop("Manifest sample column is not present in input metadata")
}

for (idx in seq_len(nrow(manifest))) {
  row <- manifest[idx, ]
  object_path <- row$object_path

  if (!file.exists(object_path)) {
    manifest_dir_path <- file.path(
      dirname(opt$manifest),
      object_path
    )

    if (!file.exists(manifest_dir_path)) {
      stop("Sample object does not exist: ", object_path)
    }

    object_path <- manifest_dir_path
  }

  sample_object <- readRDS(object_path)
  sample_values <- unique(
    as.character(sample_object[[]][[sample_column]])
  )

  if (length(sample_values) != 1 || sample_values[[1]] != row$sample) {
    stop("Sample object metadata does not match manifest for ", row$sample)
  }

  if (ncol(sample_object) != as.integer(row$cells)) {
    stop("Sample object cell count mismatch for ", row$sample)
  }

  if (nrow(sample_object) != nrow(seu)) {
    stop("Sample object feature count mismatch for ", row$sample)
  }
}

cat("Split Seurat by sample contract: VALID\n")
