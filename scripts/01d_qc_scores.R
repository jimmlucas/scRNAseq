#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
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
  ),
  make_option(
    c("--output"),
    type = "character",
    default = "seu_qc_scored.rds",
    help = "Output Seurat RDS object with QC scores"
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

get_counts_matrix <- function(seu, assay = DefaultAssay(seu)) {

  assay_obj <- seu[[assay]]

  available_layers <- Layers(assay_obj)

  message(
    "Available layers in assay ",
    assay,
    ": ",
    paste(available_layers, collapse = ", ")
  )

  # Seurat v5 can store one counts layer per sample after merge:
  # counts.wta1, counts.wta2, counts.wta3, etc.
  counts_layers <- grep("^counts", available_layers, value = TRUE)

  if (length(counts_layers) == 0) {
    stop(
      paste0(
        "No counts-like layers found in assay ",
        assay,
        ". Available layers: ",
        paste(available_layers, collapse = ", ")
      )
    )
  }

  message("Counts layers detected: ", paste(counts_layers, collapse = ", "))

  layer_mats <- lapply(counts_layers, function(layer_name) {

    message("Reading layer: ", layer_name)

    mat <- LayerData(
      object = assay_obj,
      layer = layer_name
    )

    mat <- as(mat, "dgCMatrix")

    message(
      "Layer ",
      layer_name,
      ": ",
      nrow(mat),
      " genes x ",
      ncol(mat),
      " cells"
    )

    mat
  })

  all_genes <- unique(unlist(lapply(layer_mats, rownames)))

  message("Total genes after layer union: ", length(all_genes))

  align_matrix <- function(mat, all_genes) {

    if (identical(rownames(mat), all_genes)) {
      return(mat)
    }

    aligned <- Matrix::Matrix(
      0,
      nrow = length(all_genes),
      ncol = ncol(mat),
      sparse = TRUE
    )

    rownames(aligned) <- all_genes
    colnames(aligned) <- colnames(mat)

    idx <- match(rownames(mat), all_genes)

    aligned[idx, ] <- mat

    aligned
  }

  layer_mats <- lapply(layer_mats, align_matrix, all_genes = all_genes)

  counts <- do.call(Matrix::cbind2, layer_mats)

  counts <- as(counts, "dgCMatrix")

  message(
    "Combined counts matrix: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " cells"
  )

  counts
}

match_genes <- function(gene_set, available_genes) {
  matched <- Seurat::CaseMatch(
    search = gene_set,
    match = available_genes
  )

  unique(matched)
}

calculate_sparse_log_score <- function(counts, genes, lib_size, scale_factor = 10000) {

  genes <- intersect(genes, rownames(counts))

  score <- rep(NA_real_, ncol(counts))
  names(score) <- colnames(counts)

  if (length(genes) == 0) {
    return(score)
  }

  lib_size[is.na(lib_size) | lib_size <= 0] <- 1

  sub_counts <- counts[genes, , drop = FALSE]

  norm_factor <- scale_factor / lib_size

  norm_counts <- sub_counts %*% Matrix::Diagonal(x = norm_factor)

  log_norm_counts <- log1p(norm_counts)

  score <- Matrix::colMeans(log_norm_counts)

  score <- as.numeric(score)
  names(score) <- colnames(counts)

  score
}

assign_cell_cycle_phase <- function(s_score, g2m_score) {

  phase <- rep("G1", length(s_score))

  phase[
    !is.na(s_score) &
      !is.na(g2m_score) &
      s_score > g2m_score &
      s_score > 0
  ] <- "S"

  phase[
    !is.na(s_score) &
      !is.na(g2m_score) &
      g2m_score > s_score &
      g2m_score > 0
  ] <- "G2M"

  phase[is.na(s_score) | is.na(g2m_score)] <- NA_character_

  phase
}

get_cell_cycle_genes <- function() {

  cc <- tryCatch(
    {
      get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
    },
    error = function(e) {
      NULL
    }
  )

  cc
}

# -----------------------------
# Load configuration
# -----------------------------

message_step("Loading configuration")

check_file(opt$config, "configuration file")
config <- yaml::read_yaml(opt$config)

scale_factor <- get_config_value(config, "scoring.scale_factor", 10000)

message("scale_factor: ", scale_factor)

# -----------------------------
# Load object
# -----------------------------

message_step("Loading QC-filtered Seurat object")

check_file(opt$input, "input Seurat object")

seu <- readRDS(opt$input)

message("Input cells: ", ncol(seu))
message("Input genes: ", nrow(seu))
message("Default assay: ", DefaultAssay(seu))

DefaultAssay(seu) <- "RNA"

# -----------------------------
# Extract counts
# -----------------------------

message_step("Extracting counts matrix")

counts <- get_counts_matrix(seu)

message("Counts genes: ", nrow(counts))
message("Counts cells: ", ncol(counts))

# Reorder metadata if needed to match counts columns
common_cells <- intersect(colnames(seu), colnames(counts))

if (length(common_cells) == 0) {
  stop("No shared cell names between Seurat object and counts matrix.")
}

if (length(common_cells) < ncol(seu)) {
  message(
    "Warning: counts matrix contains fewer matching cells than Seurat object. ",
    "Matching cells: ", length(common_cells),
    " / ", ncol(seu)
  )
}

# Keep only cells present in counts
seu <- subset(seu, cells = common_cells)
counts <- counts[, colnames(seu), drop = FALSE]

if ("nCount_RNA" %in% colnames(seu@meta.data)) {
  lib_size <- seu$nCount_RNA
} else {
  lib_size <- Matrix::colSums(counts)
}

names(lib_size) <- colnames(seu)

# -----------------------------
# Cell cycle score
# -----------------------------

message_step("Calculating cell cycle scores")

cc_genes <- get_cell_cycle_genes()

if (is.null(cc_genes)) {

  message("cc.genes.updated.2019 not available. Cell cycle scores will be set to NA.")

  seu$S.Score <- NA_real_
  seu$G2M.Score <- NA_real_
  seu$Phase <- NA_character_

  s_genes <- character(0)
  g2m_genes <- character(0)

} else {

  s_genes <- match_genes(
    gene_set = cc_genes$s.genes,
    available_genes = rownames(counts)
  )

  g2m_genes <- match_genes(
    gene_set = cc_genes$g2m.genes,
    available_genes = rownames(counts)
  )

  message("S genes matched: ", length(s_genes))
  message("G2M genes matched: ", length(g2m_genes))

  if (length(s_genes) < 5 || length(g2m_genes) < 5) {

    message("Too few cell cycle genes matched. Scores will be set to NA.")

    seu$S.Score <- NA_real_
    seu$G2M.Score <- NA_real_
    seu$Phase <- NA_character_

  } else {

    s_score <- calculate_sparse_log_score(
      counts = counts,
      genes = s_genes,
      lib_size = lib_size,
      scale_factor = scale_factor
    )

    g2m_score <- calculate_sparse_log_score(
      counts = counts,
      genes = g2m_genes,
      lib_size = lib_size,
      scale_factor = scale_factor
    )

    seu$S.Score <- s_score[colnames(seu)]
    seu$G2M.Score <- g2m_score[colnames(seu)]
    seu$Phase <- assign_cell_cycle_phase(seu$S.Score, seu$G2M.Score)
  }
}

gc()

# -----------------------------
# Tissue dissociation score
# -----------------------------

message_step("Calculating tissue dissociation score")

dissociation_genes <- c(
  "Fos", "Zfp36", "Fosb", "Dusp1", "Atf3", "Nr4a1", "Ppp1r15a",
  "Junb", "Egr1", "Hspa1b", "Hspa1a", "Socs3", "Klf6", "Jun",
  "Ier2", "Nfkbia", "Dnajb1", "Ier3", "Ccnl1", "Ier5", "Id1",
  "Cebpd", "Ccn1", "Irf1", "Dusp2", "Btg2", "Plaur", "Maff",
  "Klf4", "Phlda2", "Tnfaip3", "Klf10", "Ncoa7", "Hes1",
  "Zfand2a", "Ifrd1", "Ldlr", "Myc", "Mcl1", "Nr4a2", "Rhob",
  "Jund", "Pmaip1", "Nfkbiz", "Zc3h12a", "Edn1", "Sertad1",
  "Klf2", "Trib1", "Birc3", "Sox9", "Icam1", "Gadd45b",
  "Cdkn1a", "Hsph1", "Hspb1", "Emp1", "Rasd1", "Dnaja1",
  "Sat1", "Csrnp1", "Pim1", "Anxa1", "Hsp90aa1", "Areg",
  "Id2", "Tob1", "Lmna", "Gadd45a", "Dnajb4", "Rel", "Stat3",
  "Nampt", "Cebpb", "Zfp36l1", "Bag3", "Sqstm1", "Atf4"
)

dissociation_genes <- match_genes(
  gene_set = dissociation_genes,
  available_genes = rownames(counts)
)

message("Dissociation genes matched: ", length(dissociation_genes))

if (length(dissociation_genes) < 5) {

  message("Too few dissociation genes matched. Dissociation score will be set to NA.")
  seu$dissoc_score1 <- NA_real_

} else {

  dissoc_score <- calculate_sparse_log_score(
    counts = counts,
    genes = dissociation_genes,
    lib_size = lib_size,
    scale_factor = scale_factor
  )

  seu$dissoc_score1 <- dissoc_score[colnames(seu)]
}

gc()

# -----------------------------
# Summaries
# -----------------------------

message_step("Creating score summaries")

score_summary <- data.frame(
  metric = c(
    "cells",
    "genes",
    "s_genes_matched",
    "g2m_genes_matched",
    "dissociation_genes_matched",
    "median_S.Score",
    "median_G2M.Score",
    "median_dissoc_score1"
  ),
  value = c(
    ncol(seu),
    nrow(seu),
    length(s_genes),
    length(g2m_genes),
    length(dissociation_genes),
    median(seu$S.Score, na.rm = TRUE),
    median(seu$G2M.Score, na.rm = TRUE),
    median(seu$dissoc_score1, na.rm = TRUE)
  )
)

readr::write_csv(score_summary, "01d_qc_scores_summary.csv")

phase_summary <- seu@meta.data %>%
  dplyr::count(Phase, name = "cells") %>%
  dplyr::arrange(Phase)

readr::write_csv(phase_summary, "01d_cell_cycle_phase_summary.csv")

if ("sample" %in% colnames(seu@meta.data)) {

  if (!"Group" %in% colnames(seu@meta.data)) {
    seu$Group <- NA_character_
  }

  score_by_sample <- seu@meta.data %>%
    dplyr::group_by(sample, Group) %>%
    dplyr::summarise(
      cells = dplyr::n(),
      median_S.Score = median(S.Score, na.rm = TRUE),
      median_G2M.Score = median(G2M.Score, na.rm = TRUE),
      median_dissoc_score1 = median(dissoc_score1, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(score_by_sample, "01d_qc_scores_by_sample.csv")
}

score_metadata <- seu@meta.data
score_metadata$cell_id <- rownames(score_metadata)

readr::write_csv(score_metadata, "01d_qc_scores_metadata.csv")

# -----------------------------
# Save object
# -----------------------------

message_step("Saving outputs")

saveRDS(seu, opt$output)

sink("01d_sessionInfo.txt")
print(sessionInfo())
sink()

message("Scored object saved in: ", opt$output)
message("Score summary saved in: 01d_qc_scores_summary.csv")
message("Phase summary saved in: 01d_cell_cycle_phase_summary.csv")
message("Score metadata saved in: 01d_qc_scores_metadata.csv")
message("Session info saved in: 01d_sessionInfo.txt")

message_step("01d_qc_scores finished successfully")