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
    c("--input"),
    type = "character",
    default = NULL,
    help = "Input Seurat RDS object from 00_load_data"
  ),
  make_option(
    c("--config"),
    type = "character",
    default = NULL,
    help = "YAML configuration file"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = "seu_qc_filtered.rds",
    help = "Output filtered Seurat RDS object"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Missing --input")
}

if (is.null(opt$config)) {
  stop("Missing --config")
}

# -----------------------------
# Helper functions
# -----------------------------

message_step <- function(x) {
  message("\n========== ", x, " ==========")
}

check_file <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(paste0("File does not exist - ", label, ": ", path))
  }
}

get_config_value <- function(config, path, default = NULL) {
  parts <- strsplit(path, "\\.")[[1]]
  value <- config

  for (p in parts) {
    if (is.null(value[[p]])) {
      return(default)
    }
    value <- value[[p]]
  }

  value
}

safe_percentage_feature_set <- function(seu, pattern, colname) {
  genes <- grep(pattern, rownames(seu), value = TRUE)

  if (length(genes) == 0) {
    message("No genes found for ", colname, " using pattern: ", pattern)
    seu[[colname]] <- 0
  } else {
    message(colname, " genes detected: ", length(genes))
    seu[[colname]] <- PercentageFeatureSet(
      object = seu,
      features = genes
    )
  }

  seu
}

# -----------------------------
# Load configuration
# -----------------------------

message_step("Loading configuration")

check_file(opt$config, "configuration file")
config <- yaml::read_yaml(opt$config)

organism <- get_config_value(config, "project.organism", "mouse")

# QC thresholds adapted from the 01a QC filtering step
min_features <- get_config_value(config, "qc.min_features", 1000)
max_features <- get_config_value(config, "qc.max_features", 100000)
min_counts   <- get_config_value(config, "qc.min_counts", 0)
max_counts   <- get_config_value(config, "qc.max_counts", 40000)
max_mito     <- get_config_value(config, "qc.max_percent_mito", 20)
max_ribo     <- get_config_value(config, "qc.max_percent_ribo", 15)
max_globin   <- get_config_value(config, "qc.max_percent_globin", 0.1)

message("organism: ", organism)
message("min_features: ", min_features)
message("max_features: ", max_features)
message("min_counts: ", min_counts)
message("max_counts: ", max_counts)
message("max_percent_mito: ", max_mito)
message("max_percent_ribo: ", max_ribo)
message("max_percent_globin: ", max_globin)

# -----------------------------
# Load Seurat object
# -----------------------------

message_step("Loading Seurat object")

check_file(opt$input, "input Seurat object")

seu <- readRDS(opt$input)

message("Input cells: ", ncol(seu))
message("Input genes: ", nrow(seu))
message("Default assay: ", DefaultAssay(seu))

DefaultAssay(seu) <- "RNA"

# -----------------------------
# Ensure metadata columns
# -----------------------------

message_step("Checking metadata")

if (!"sample" %in% colnames(seu@meta.data)) {
  if ("ShortName" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$ShortName
  } else if ("sample_id" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$sample_id
  } else if ("orig.ident" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$orig.ident
  } else {
    seu$sample <- "sample"
  }
}

if (!"Group" %in% colnames(seu@meta.data)) {
  seu$Group <- NA_character_
}

# -----------------------------
# Calculate QC metrics
# -----------------------------

message_step("Calculating QC metrics")

if (tolower(organism) %in% c("mouse", "mus_musculus", "mmusculus")) {
  mito_pattern   <- "^mt-"
  ribo_pattern   <- "^Rpl|^Rps"
  globin_pattern <- "^Hba|^Hbb|^Hb[ab]"
} else {
  mito_pattern   <- "^MT-"
  ribo_pattern   <- "^RPL|^RPS"
  globin_pattern <- "^HBA|^HBB|^HB[AB]"
}

seu <- safe_percentage_feature_set(seu, mito_pattern, "percent.mito")
seu <- safe_percentage_feature_set(seu, ribo_pattern, "percent.ribo")
seu <- safe_percentage_feature_set(seu, globin_pattern, "percent.globin")

seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
seu$cell_id <- colnames(seu)

# -----------------------------
# Summary before filtering
# -----------------------------

message_step("Summary before filtering")

cells_before <- ncol(seu)
genes_before <- nrow(seu)

summary_before <- data.frame(
  metric = c(
    "cells_before_filtering",
    "genes_before_filtering",
    "median_nCount_RNA_before",
    "median_nFeature_RNA_before",
    "median_percent_mito_before",
    "median_percent_ribo_before",
    "median_percent_globin_before",
    "median_log10GenesPerUMI_before"
  ),
  value = c(
    cells_before,
    genes_before,
    median(seu$nCount_RNA, na.rm = TRUE),
    median(seu$nFeature_RNA, na.rm = TRUE),
    median(seu$percent.mito, na.rm = TRUE),
    median(seu$percent.ribo, na.rm = TRUE),
    median(seu$percent.globin, na.rm = TRUE),
    median(seu$log10GenesPerUMI, na.rm = TRUE)
  )
)

cells_by_sample_before <- seu@meta.data %>%
  dplyr::group_by(sample, Group) %>%
  dplyr::summarise(
    cells_before = dplyr::n(),
    median_nCount_RNA_before = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA_before = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mito_before = median(percent.mito, na.rm = TRUE),
    median_percent_ribo_before = median(percent.ribo, na.rm = TRUE),
    median_percent_globin_before = median(percent.globin, na.rm = TRUE),
    median_log10GenesPerUMI_before = median(log10GenesPerUMI, na.rm = TRUE),
    .groups = "drop"
  )

message("Cells before filtering: ", cells_before)
message("Genes before filtering: ", genes_before)

# -----------------------------
# Apply QC filters
# -----------------------------

message_step("Applying QC filters")

seu <- subset(
  x = seu,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA <= max_features &
    nCount_RNA >= min_counts &
    nCount_RNA < max_counts &
    percent.mito < max_mito &
    percent.ribo < max_ribo &
    percent.globin < max_globin
)

gc()

cells_after <- ncol(seu)
genes_after <- nrow(seu)

message("Cells after filtering: ", cells_after)
message("Genes after filtering: ", genes_after)
message("Cells removed: ", cells_before - cells_after)

if (cells_after == 0) {
  stop("No cells remained after QC filtering. Relax QC thresholds.")
}

# -----------------------------
# Summary after filtering
# -----------------------------

message_step("Summary after filtering")

summary_after <- data.frame(
  metric = c(
    "cells_after_filtering",
    "genes_after_filtering",
    "cells_removed",
    "median_nCount_RNA_after",
    "median_nFeature_RNA_after",
    "median_percent_mito_after",
    "median_percent_ribo_after",
    "median_percent_globin_after",
    "median_log10GenesPerUMI_after"
  ),
  value = c(
    cells_after,
    genes_after,
    cells_before - cells_after,
    median(seu$nCount_RNA, na.rm = TRUE),
    median(seu$nFeature_RNA, na.rm = TRUE),
    median(seu$percent.mito, na.rm = TRUE),
    median(seu$percent.ribo, na.rm = TRUE),
    median(seu$percent.globin, na.rm = TRUE),
    median(seu$log10GenesPerUMI, na.rm = TRUE)
  )
)

cells_by_sample_after <- seu@meta.data %>%
  dplyr::group_by(sample, Group) %>%
  dplyr::summarise(
    cells_after = dplyr::n(),
    median_nCount_RNA_after = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA_after = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mito_after = median(percent.mito, na.rm = TRUE),
    median_percent_ribo_after = median(percent.ribo, na.rm = TRUE),
    median_percent_globin_after = median(percent.globin, na.rm = TRUE),
    median_log10GenesPerUMI_after = median(log10GenesPerUMI, na.rm = TRUE),
    .groups = "drop"
  )

qc_by_sample <- dplyr::full_join(
  cells_by_sample_before,
  cells_by_sample_after,
  by = c("sample", "Group")
) %>%
  dplyr::mutate(
    cells_removed = cells_before - cells_after
  )

qc_summary <- dplyr::bind_rows(summary_before, summary_after)

# -----------------------------
# Placeholder score columns
# -----------------------------
# Scores are intentionally not calculated here.
# Cell cycle and dissociation scores should be handled in a later module.

message_step("Adding placeholder score columns")

seu$S.Score <- NA_real_
seu$G2M.Score <- NA_real_
seu$Phase <- NA_character_
seu$dissoc_score1 <- NA_real_

gc()

# -----------------------------
# Save outputs
# -----------------------------

message_step("Saving outputs")

saveRDS(seu, opt$output)

# Required by Nextflow
readr::write_csv(qc_summary, "01_qc_summary.csv")
sink("01_sessionInfo.txt")
print(sessionInfo())
sink()

# Optional diagnostic outputs
readr::write_csv(qc_summary, "01_qc_filter_summary.csv")
readr::write_csv(qc_by_sample, "01_qc_filter_by_sample.csv")

cell_metadata <- seu@meta.data
cell_metadata$cell_id <- rownames(cell_metadata)
readr::write_csv(cell_metadata, "01_qc_cell_metadata.csv")

message("Filtered object saved in: ", opt$output)
message("QC summary saved in: 01_qc_summary.csv")
message("QC summary copy saved in: 01_qc_filter_summary.csv")
message("QC by sample saved in: 01_qc_filter_by_sample.csv")
message("Cell metadata saved in: 01_qc_cell_metadata.csv")
message("Session info saved in: 01_sessionInfo.txt")

message_step("01_qc_filter finished successfully")
