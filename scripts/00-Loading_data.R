#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -----------------------------
# Arguments
# -----------------------------

option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file"
  ),
  make_option(
    c("--project_dir"),
    type = "character",
    default = NULL,
    help = "Path to the Nextflow project directory"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config)) {
  stop("You must specify a configuration file with --config")
}

if (is.null(opt$project_dir)) {
  stop("You must specify the project directory with --project_dir")
}

config <- yaml::read_yaml(opt$config)

# -----------------------------
# Auxiliary functions
# -----------------------------

check_file <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(paste0("File does not exist - ", label, ": ", path))
  }
}

check_dir <- function(path, label = "directory") {
  if (!dir.exists(path)) {
    stop(paste0("Directory does not exist - ", label, ": ", path))
  }
}

message_step <- function(x) {
  message("\n========== ", x, " ==========")
}

find_seurat_rds <- function(sample_dir, sample_id) {
  rds_files <- list.files(
    sample_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(rds_files) == 0) {
    stop("No .rds file found for sample ", sample_id, " in: ", sample_dir)
  }

  if (length(rds_files) > 1) {
    stop(
      paste0(
        "Multiple .rds files found for sample ", sample_id, " in: ", sample_dir, "\n",
        "Files found:\n",
        paste(rds_files, collapse = "\n"),
        "\nKeep only one .rds file per sample directory."
      )
    )
  }

  rds_files[1]
}

get_counts_matrix <- function(obj) {
  assay_name <- DefaultAssay(obj)

  counts <- tryCatch(
    {
      GetAssayData(obj, assay = assay_name, layer = "counts")
    },
    error = function(e) {
      GetAssayData(obj, assay = assay_name, slot = "counts")
    }
  )

  counts
}

make_minimal_seurat_object <- function(obj, sample_id, row_meta) {
  counts <- get_counts_matrix(obj)

  minimal_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id
  )

  minimal_obj$sample_id <- sample_id
  minimal_obj$orig.ident <- sample_id

  for (col in colnames(row_meta)) {
    value <- row_meta[[col]][1]
    minimal_obj[[col]] <- value
  }

  minimal_obj <- RenameCells(minimal_obj, add.cell.id = sample_id)

  minimal_obj
}

# -----------------------------
# Load configuration
# -----------------------------

project_dir <- normalizePath(opt$project_dir, mustWork = TRUE)

project_name <- config$project$name

targets_file <- file.path(project_dir, config$input$targets)
count_dir    <- file.path(project_dir, config$input$count_dir)

sample_column <- config$project$sample_column
short_column  <- config$project$shortname_column
group_column  <- config$project$group_column

# Nextflow expects these files in the current work directory
output_rds    <- config$output$merged_object
summary_file  <- "00_load_data_summary.csv"
session_file  <- "00_sessionInfo.txt"

message_step("Project paths")

message("Project directory: ", project_dir)
message("Targets file: ", targets_file)
message("Count directory: ", count_dir)
message("Output RDS: ", output_rds)

# -----------------------------
# Validate inputs
# -----------------------------

message_step("Validating inputs")

check_file(targets_file, "metadata file")
check_dir(count_dir, "counts directory")

# -----------------------------
# Read metadata table
# -----------------------------

message_step("Reading samples table")

targets <- readr::read_csv(targets_file, show_col_types = FALSE)

required_cols <- c(sample_column, short_column)

if (!is.null(group_column)) {
  required_cols <- c(required_cols, group_column)
}

missing_cols <- setdiff(required_cols, colnames(targets))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns in targets: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (anyDuplicated(targets[[short_column]]) > 0) {
  stop("There are duplicate values in the short identifier column: ", short_column)
}

print(targets)

# -----------------------------
# Resolve Seurat object paths
# -----------------------------

message_step("Resolving Seurat objects")

sample_ids <- targets[[short_column]]
sample_dirs <- file.path(count_dir, sample_ids)
names(sample_dirs) <- sample_ids

for (nm in names(sample_dirs)) {
  check_dir(sample_dirs[[nm]], paste0("sample directory for ", nm))
}

seurat_paths <- vapply(
  names(sample_dirs),
  function(nm) {
    find_seurat_rds(sample_dirs[[nm]], nm)
  },
  character(1)
)

names(seurat_paths) <- sample_ids

message("Resolved Seurat RDS files:")
print(seurat_paths)

# -----------------------------
# Load, rebuild minimal objects and merge incrementally
# -----------------------------

message_step("Loading minimal Seurat objects and merging incrementally")

summary_list <- list()
seu <- NULL

for (nm in names(seurat_paths)) {

  message_step(paste0("Processing sample: ", nm))

  row_meta <- targets[targets[[short_column]] == nm, , drop = FALSE]

  if (nrow(row_meta) != 1) {
    stop("Could not find a unique metadata row for: ", nm)
  }

  message("Reading: ", seurat_paths[[nm]])

  raw_obj <- readRDS(seurat_paths[[nm]])

  message("Original cells: ", ncol(raw_obj))
  message("Original genes: ", nrow(raw_obj))
  message("Original default assay: ", DefaultAssay(raw_obj))

  obj <- make_minimal_seurat_object(
    obj = raw_obj,
    sample_id = nm,
    row_meta = row_meta
  )

  rm(raw_obj)
  gc()

  message("Minimal object cells: ", ncol(obj))
  message("Minimal object genes: ", nrow(obj))

  summary_list[[nm]] <- data.frame(
    sample_id = nm,
    cells = ncol(obj),
    genes = nrow(obj),
    seurat_file = seurat_paths[[nm]],
    stringsAsFactors = FALSE
  )

  if (is.null(seu)) {

    seu <- obj

  } else {

    message("Merging sample into combined object: ", nm)

    seu <- merge(
      x = seu,
      y = obj,
      project = project_name,
      merge.data = FALSE,
      merge.dr = FALSE
    )
  }

  rm(obj)
  gc()

  message("Current merged cells: ", ncol(seu))
  message("Current merged genes: ", nrow(seu))
}

if (is.null(seu)) {
  stop("No Seurat objects were loaded.")
}

message_step("Final merged object")

message("Merged object cells: ", ncol(seu))
message("Merged object genes: ", nrow(seu))
message("Default assay: ", DefaultAssay(seu))

# -----------------------------
# Save merged object
# -----------------------------

message_step("Saving combined object")

saveRDS(seu, output_rds)

message("Object saved in: ", output_rds)

# -----------------------------
# Save summary
# -----------------------------

summary_df <- dplyr::bind_rows(summary_list)

readr::write_csv(summary_df, summary_file)

message("Summary saved in: ", summary_file)

# -----------------------------
# Session info
# -----------------------------

sink(session_file)
print(sessionInfo())
sink()

message("Session info saved in: ", session_file)

message_step("00_load_data finished successfully")
