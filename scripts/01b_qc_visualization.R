#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# -----------------------------
# Arguments
# -----------------------------

option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "Input QC-filtered Seurat RDS object"
  ),
  make_option(
    c("--config"),
    type = "character",
    default = NULL,
    help = "YAML configuration file"
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

safe_save_plot <- function(plot, filename, width = 8, height = 5, dpi = 150) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  message("Plot saved: ", filename)
}

# -----------------------------
# Load config
# -----------------------------

message_step("Loading configuration")

check_file(opt$config, "configuration file")
config <- yaml::read_yaml(opt$config)

project_name <- get_config_value(config, "project.name", "scRNAseq_project")

message("Project: ", project_name)

# -----------------------------
# Load object
# -----------------------------

message_step("Loading QC-filtered Seurat object")

check_file(opt$input, "input Seurat object")

seu <- readRDS(opt$input)

message("Input cells: ", ncol(seu))
message("Input genes: ", nrow(seu))
message("Default assay: ", DefaultAssay(seu))

# -----------------------------
# Validate metadata
# -----------------------------

message_step("Checking metadata columns")

required_cols <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mito",
  "percent.ribo",
  "percent.globin",
  "log10GenesPerUMI"
)

missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required QC metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

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
# Export metadata
# -----------------------------

message_step("Exporting QC metadata")

qc_meta <- seu@meta.data
qc_meta$cell_id <- rownames(qc_meta)

readr::write_csv(qc_meta, "01b_qc_cell_metadata.csv")

qc_summary <- data.frame(
  metric = c(
    "cells",
    "genes",
    "median_nCount_RNA",
    "median_nFeature_RNA",
    "median_percent_mito",
    "median_percent_ribo",
    "median_percent_globin",
    "median_log10GenesPerUMI"
  ),
  value = c(
    ncol(seu),
    nrow(seu),
    median(seu$nCount_RNA, na.rm = TRUE),
    median(seu$nFeature_RNA, na.rm = TRUE),
    median(seu$percent.mito, na.rm = TRUE),
    median(seu$percent.ribo, na.rm = TRUE),
    median(seu$percent.globin, na.rm = TRUE),
    median(seu$log10GenesPerUMI, na.rm = TRUE)
  )
)

readr::write_csv(qc_summary, "01b_qc_summary.csv")

qc_by_sample <- qc_meta %>%
  dplyr::group_by(sample, Group) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mito = median(percent.mito, na.rm = TRUE),
    median_percent_ribo = median(percent.ribo, na.rm = TRUE),
    median_percent_globin = median(percent.globin, na.rm = TRUE),
    median_log10GenesPerUMI = median(log10GenesPerUMI, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(qc_by_sample, "01b_qc_by_sample.csv")

qc_by_group <- qc_meta %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mito = median(percent.mito, na.rm = TRUE),
    median_percent_ribo = median(percent.ribo, na.rm = TRUE),
    median_percent_globin = median(percent.globin, na.rm = TRUE),
    median_log10GenesPerUMI = median(log10GenesPerUMI, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(qc_by_group, "01b_qc_by_group.csv")

# -----------------------------
# QC plots by sample
# -----------------------------

message_step("Generating QC plots by sample")

p_ncount <- ggplot(qc_meta, aes(x = sample, y = nCount_RNA)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "UMI counts per cell after QC",
    x = "Sample",
    y = "nCount_RNA"
  )

safe_save_plot(p_ncount, "01b_qc_nCount_RNA_by_sample.png")

p_nfeature <- ggplot(qc_meta, aes(x = sample, y = nFeature_RNA)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Detected genes per cell after QC",
    x = "Sample",
    y = "nFeature_RNA"
  )

safe_save_plot(p_nfeature, "01b_qc_nFeature_RNA_by_sample.png")

p_mito <- ggplot(qc_meta, aes(x = sample, y = percent.mito)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mitochondrial percentage after QC",
    x = "Sample",
    y = "percent.mito"
  )

safe_save_plot(p_mito, "01b_qc_percent_mito_by_sample.png")

p_ribo <- ggplot(qc_meta, aes(x = sample, y = percent.ribo)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Ribosomal percentage after QC",
    x = "Sample",
    y = "percent.ribo"
  )

safe_save_plot(p_ribo, "01b_qc_percent_ribo_by_sample.png")

p_globin <- ggplot(qc_meta, aes(x = sample, y = percent.globin)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Globin percentage after QC",
    x = "Sample",
    y = "percent.globin"
  )

safe_save_plot(p_globin, "01b_qc_percent_globin_by_sample.png")

p_complexity <- ggplot(qc_meta, aes(x = sample, y = log10GenesPerUMI)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Library complexity after QC",
    x = "Sample",
    y = "log10GenesPerUMI"
  )

safe_save_plot(p_complexity, "01b_qc_log10GenesPerUMI_by_sample.png")

# -----------------------------
# QC scatter plots
# -----------------------------

message_step("Generating QC scatter plots")

p_scatter_complexity <- ggplot(
  qc_meta,
  aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mito)
) +
  geom_point(size = 0.35, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Cell complexity after QC",
    x = "nCount_RNA",
    y = "nFeature_RNA",
    colour = "percent.mito"
  )

safe_save_plot(
  p_scatter_complexity,
  "01b_qc_scatter_nCount_nFeature_percent_mito.png",
  width = 7,
  height = 5
)

p_scatter_mito <- ggplot(
  qc_meta,
  aes(x = nCount_RNA, y = percent.mito)
) +
  geom_point(size = 0.35, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Mitochondrial percentage vs UMI counts",
    x = "nCount_RNA",
    y = "percent.mito"
  )

safe_save_plot(
  p_scatter_mito,
  "01b_qc_scatter_nCount_percent_mito.png",
  width = 7,
  height = 5
)

p_scatter_ribo <- ggplot(
  qc_meta,
  aes(x = nCount_RNA, y = percent.ribo)
) +
  geom_point(size = 0.35, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Ribosomal percentage vs UMI counts",
    x = "nCount_RNA",
    y = "percent.ribo"
  )

safe_save_plot(
  p_scatter_ribo,
  "01b_qc_scatter_nCount_percent_ribo.png",
  width = 7,
  height = 5
)

# -----------------------------
# Group-level plots
# -----------------------------

message_step("Generating QC plots by group")

if (length(unique(na.omit(qc_meta$Group))) > 1) {

  p_group_ncount <- ggplot(qc_meta, aes(x = Group, y = nCount_RNA)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.2) +
    theme_bw() +
    labs(
      title = "UMI counts per group after QC",
      x = "Group",
      y = "nCount_RNA"
    )

  safe_save_plot(p_group_ncount, "01b_qc_nCount_RNA_by_group.png")

  p_group_nfeature <- ggplot(qc_meta, aes(x = Group, y = nFeature_RNA)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.2) +
    theme_bw() +
    labs(
      title = "Detected genes per group after QC",
      x = "Group",
      y = "nFeature_RNA"
    )

  safe_save_plot(p_group_nfeature, "01b_qc_nFeature_RNA_by_group.png")

  p_group_mito <- ggplot(qc_meta, aes(x = Group, y = percent.mito)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.2) +
    theme_bw() +
    labs(
      title = "Mitochondrial percentage per group after QC",
      x = "Group",
      y = "percent.mito"
    )

  safe_save_plot(p_group_mito, "01b_qc_percent_mito_by_group.png")
}

# -----------------------------
# Save session info
# -----------------------------

message_step("Saving session info")

sink("01b_sessionInfo.txt")
print(sessionInfo())
sink()

message("QC visualization summary saved in: 01b_qc_summary.csv")
message("QC by sample saved in: 01b_qc_by_sample.csv")
message("QC by group saved in: 01b_qc_by_group.csv")
message("QC metadata saved in: 01b_qc_cell_metadata.csv")
message("Session info saved in: 01b_sessionInfo.txt")

message_step("01b_qc_visualization finished successfully")