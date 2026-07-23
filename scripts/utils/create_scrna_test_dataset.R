#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

find_single_rds <- function(sample_dir) {
  rds_files <- list.files(
    sample_dir,
    pattern = "_Seurat\\.rds$|\\.rds$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  if (length(rds_files) == 0) {
    stop("No Seurat .rds file found in: ", sample_dir)
  }

  if (length(rds_files) > 1) {
    stop(
      "Multiple .rds files found in: ",
      sample_dir,
      "\n",
      paste(rds_files, collapse = "\n")
    )
  }

  rds_files[1]
}

make_small_seurat <- function(sample_dir, output_rds, sample_name, short_name, group, n_cells = 200) {
  input_rds <- find_single_rds(sample_dir)

  message("Reading: ", input_rds)

  obj <- readRDS(input_rds)

  message("Original dimensions: ", nrow(obj), " genes x ", ncol(obj), " cells")

  if ("RNA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "RNA"
  }

  cells <- colnames(obj)
  n_take <- min(n_cells, length(cells))

  selected_cells <- sample(cells, n_take)

  obj <- subset(obj, cells = selected_cells)

  obj$SampleName <- sample_name
  obj$ShortName <- short_name
  obj$Group <- group
  obj$sample_id <- short_name
  obj$sample <- short_name
  obj$orig.ident <- short_name

  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()
  obj@commands <- list()

  message("Test dimensions: ", nrow(obj), " genes x ", ncol(obj), " cells")

  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

  saveRDS(obj, output_rds)

  message("Saved: ", output_rds)

  invisible(NULL)
}

make_small_seurat(
  sample_dir = "data/count/wta1",
  output_rds = "data_test/count/wta1/wta1_test_Seurat.rds",
  sample_name = "wta-1",
  short_name = "wta1",
  group = "CTL",
  n_cells = 200
)

make_small_seurat(
  sample_dir = "data/count/wta5",
  output_rds = "data_test/count/wta5/wta5_test_Seurat.rds",
  sample_name = "wta-5",
  short_name = "wta5",
  group = "HF",
  n_cells = 200
)

targets <- data.frame(
  SampleName = c("wta-1", "wta-5"),
  ShortName = c("wta1", "wta5"),
  Group = c("CTL", "HF")
)

dir.create("data_test", recursive = TRUE, showWarnings = FALSE)

write.csv(
  targets,
  "data_test/targets_test.csv",
  row.names = FALSE,
  quote = FALSE
)

message("Saved: data_test/targets_test.csv")
message("scRNA-seq test dataset created successfully.")
