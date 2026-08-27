#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

# ==============================================================================
# 00-Loading_data.R
#
# Canonical loader for scRNAseq.
#
# Supported inputs:
#   - bd_rhapsody : primary BD Rhapsody matrices + vendor metadata
#   - seurat_rds  : legacy/pre-existing Seurat RDS objects
#
# RESPONSIBILITIES
#   ✓ load counts
#   ✓ preserve vendor metadata
#   ✓ preserve Sample Tag calls
#   ✓ preserve Multiplet / Undetermined
#   ✓ preserve optional AB/ADT assay
#   ✓ add physical capture identity
#   ✓ map capture + Sample Tag -> biological sample metadata
#   ✓ generate globally unique cell IDs
#   ✓ merge incrementally
#   ✓ write provenance summaries
#
# THIS SCRIPT MUST NOT
#   ✗ remove cells
#   ✗ remove multiplets
#   ✗ remove undetermined cells
#   ✗ perform QC filtering
#   ✗ normalize
#   ✗ run PCA/UMAP/clustering
#   ✗ perform transcriptomic doublet detection
# ==============================================================================


# ==============================================================================
# CLI
# ==============================================================================

option_list <- list(

  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to YAML configuration file"
  ),

  make_option(
    c("--project_dir"),
    type = "character",
    default = NULL,
    help = "Nextflow project directory"
  )
)

opt <- parse_args(
  OptionParser(
    option_list = option_list
  )
)

if (is.null(opt$config)) {
  stop("You must specify --config")
}

if (is.null(opt$project_dir)) {
  stop("You must specify --project_dir")
}


# ==============================================================================
# GENERAL HELPERS
# ==============================================================================

message_step <- function(x) {
  message(
    "\n========== ",
    x,
    " =========="
  )
}


`%||%` <- function(x, y) {

  if (
    is.null(x) ||
    length(x) == 0
  ) {
    return(y)
  }

  x
}


is_missing_value <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0
  ) {
    return(TRUE)
  }

  x <- as.character(x[1])

  is.na(x) ||
    trimws(x) == "" ||
    tolower(trimws(x)) == "na"
}


check_file <- function(
  path,
  label = "file"
) {

  if (!file.exists(path)) {

    stop(
      "File does not exist - ",
      label,
      ": ",
      path
    )
  }
}


check_dir <- function(
  path,
  label = "directory"
) {

  if (!dir.exists(path)) {

    stop(
      "Directory does not exist - ",
      label,
      ": ",
      path
    )
  }
}


resolve_path <- function(
  project_dir,
  path
) {

  if (is_missing_value(path)) {
    return(NA_character_)
  }

  path <- as.character(path[1])

  if (grepl("^/", path)) {
    return(path)
  }

  file.path(
    project_dir,
    path
  )
}


safe_character <- function(x) {

  if (is.factor(x)) {
    return(as.character(x))
  }

  x
}


# ==============================================================================
# CONFIGURATION
# ==============================================================================

config <- yaml::read_yaml(
  opt$config
)

project_dir <- normalizePath(
  opt$project_dir,
  mustWork = TRUE
)

project_name <-
  config$project$name %||%
  "scRNAseq"

input_type <-
  config$input$type %||%
  "seurat_rds"

targets_file <- resolve_path(
  project_dir,
  config$input$targets
)

#
# Important:
# modules/load.nf currently expects exactly "seu.rds".
#

output_rds <- "seu.rds"

summary_file <-
  "00_load_data_summary.csv"

capture_summary_file <-
  "00_capture_summary.csv"

assignment_summary_file <-
  "00_rhapsody_assignment_summary.csv"

unmapped_assignments_file <-
  "00_unmapped_assignments.csv"

session_file <-
  "00_sessionInfo.txt"


# ==============================================================================
# TARGETS
# ==============================================================================

message_step(
  "Reading project configuration"
)

message(
  "Project: ",
  project_name
)

message(
  "Input type: ",
  input_type
)

message(
  "Targets: ",
  targets_file
)

check_file(
  targets_file,
  "targets metadata"
)

targets <- readr::read_csv(
  targets_file,
  show_col_types = FALSE
)

sample_column <-
  config$project$sample_column

shortname_column <-
  config$project$shortname_column

group_column <-
  config$project$group_column

tag_column <-
  config$project$tag_column

capture_column <-
  config$project$capture_column


required_target_columns <- c(
  sample_column,
  shortname_column
)

required_target_columns <-
  required_target_columns[
    !vapply(
      required_target_columns,
      is.null,
      logical(1)
    )
  ]

missing_target_columns <- setdiff(
  required_target_columns,
  colnames(targets)
)

if (length(missing_target_columns) > 0) {

  stop(
    "Missing required columns in targets: ",
    paste(
      missing_target_columns,
      collapse = ", "
    )
  )
}


if (input_type == "bd_rhapsody") {

  if (is.null(tag_column)) {
    stop(
      "project.tag_column must be defined for bd_rhapsody input."
    )
  }

  if (is.null(capture_column)) {
    stop(
      "project.capture_column must be defined for bd_rhapsody input."
    )
  }

  missing_rhapsody_targets <- setdiff(
    c(
      tag_column,
      capture_column
    ),
    colnames(targets)
  )

  if (
    length(
      missing_rhapsody_targets
    ) > 0
  ) {

    stop(
      "Missing Rhapsody mapping columns in targets: ",
      paste(
        missing_rhapsody_targets,
        collapse = ", "
      )
    )
  }
}


message(
  "Targets rows: ",
  nrow(targets)
)

message(
  "Targets columns: ",
  paste(
    colnames(targets),
    collapse = ", "
  )
)


# ==============================================================================
# MATRIX HELPERS
# ==============================================================================

read_primary_matrix <- function(path) {

  message(
    "Reading primary expression matrix: ",
    path
  )

  if (dir.exists(path)) {

    return(
      Seurat::Read10X(
        data.dir = path,
        unique.features = TRUE
      )
    )
  }

  if (
    file.exists(path) &&
    grepl(
      "\\.(h5|hdf5)$",
      path,
      ignore.case = TRUE
    )
  ) {

    return(
      Seurat::Read10X_h5(
        filename = path,
        unique.features = TRUE
      )
    )
  }

  stop(
    paste0(
      "Unsupported matrix input: ",
      path,
      "\n",
      "Expected either a MEX directory or an H5 file."
    )
  )
}


extract_gene_expression <- function(x) {

  if (
    inherits(x, "Matrix") ||
    inherits(x, "matrix")
  ) {
    return(x)
  }

  if (!is.list(x)) {

    stop(
      "Unsupported matrix object returned by input reader."
    )
  }

  candidates <- c(
    "Gene Expression",
    "GeneExpression",
    "RNA",
    "mRNA"
  )

  hit <- intersect(
    candidates,
    names(x)
  )

  if (length(hit) == 0) {

    stop(
      paste0(
        "Could not identify Gene Expression matrix.\n",
        "Available matrices: ",
        paste(
          names(x),
          collapse = ", "
        )
      )
    )
  }

  x[[hit[1]]]
}


extract_antibody_capture <- function(x) {

  if (!is.list(x)) {
    return(NULL)
  }

  candidates <- c(
    "Antibody Capture",
    "AntibodyCapture",
    "AbSeq",
    "AB",
    "ADT"
  )

  hit <- intersect(
    candidates,
    names(x)
  )

  if (length(hit) == 0) {
    return(NULL)
  }

  x[[hit[1]]]
}


# ==============================================================================
# BD RHAPSODY METADATA READERS
# ==============================================================================

read_rhapsody_sample_tags <- function(
  path,
  skip_lines = 7
) {

  check_file(
    path,
    "Rhapsody Sample Tag calls"
  )

  message(
    "Reading Sample Tag calls: ",
    path
  )

  x <- tryCatch(
    {

      read.csv(
        path,
        skip = skip_lines,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

    },
    error = function(e) {

      message(
        "skip=",
        skip_lines,
        " failed; trying regular CSV."
      )

      NULL
    }
  )

  if (
    is.null(x) ||
    nrow(x) == 0 ||
    ncol(x) == 0
  ) {

    x <- read.csv(
      path,
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (anyDuplicated(rownames(x))) {

    stop(
      "Duplicated barcodes in Sample Tag metadata: ",
      path
    )
  }

  x
}


read_optional_metadata <- function(path) {

  if (is_missing_value(path)) {
    return(NULL)
  }

  check_file(
    path,
    "optional cell metadata"
  )

  message(
    "Reading optional cell metadata: ",
    path
  )

  x <- read.csv(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (anyDuplicated(rownames(x))) {

    stop(
      "Duplicated barcodes in optional metadata: ",
      path
    )
  }

  x
}


# ==============================================================================
# BARCODE RECONCILIATION
# ==============================================================================

align_metadata_to_cells <- function(
  metadata,
  cells,
  label,
  barcode_prefix = NA_character_
) {

  if (is.null(metadata)) {
    return(NULL)
  }

  metadata_ids <- rownames(metadata)

  if (is.null(metadata_ids)) {

    stop(
      label,
      " does not contain barcode row names."
    )
  }


  # ---------------------------------------------------------------------------
  # 1. Exact correspondence
  # ---------------------------------------------------------------------------

  if (
    length(metadata_ids) == length(cells) &&
    setequal(
      metadata_ids,
      cells
    )
  ) {

    metadata <- metadata[
      cells,
      ,
      drop = FALSE
    ]

    return(metadata)
  }


  # ---------------------------------------------------------------------------
  # 2. Explicit configured prefix
  #
  # Historical BD exports may use:
  #
  # metadata:
  #   AAAC...
  #
  # expression:
  #   Run1_AAAC...
  #
  # This transformation is only performed when barcode_prefix is explicitly
  # supplied in the capture manifest.
  # ---------------------------------------------------------------------------

  if (!is_missing_value(barcode_prefix)) {

    prefixed_ids <- paste0(
      barcode_prefix,
      "_",
      metadata_ids
    )

    if (
      length(prefixed_ids) == length(cells) &&
      setequal(
        prefixed_ids,
        cells
      )
    ) {

      rownames(metadata) <- prefixed_ids

      metadata <- metadata[
        cells,
        ,
        drop = FALSE
      ]

      message(
        "Barcode reconciliation for ",
        label,
        ": applied explicit prefix '",
        barcode_prefix,
        "_'."
      )

      return(metadata)
    }
  }


  # ---------------------------------------------------------------------------
  # Failure: never guess barcode relationships.
  # ---------------------------------------------------------------------------

  missing_in_metadata <- setdiff(
    cells,
    metadata_ids
  )

  metadata_not_in_matrix <- setdiff(
    metadata_ids,
    cells
  )

  stop(
    paste0(
      "\nBarcode mismatch for ",
      label,
      ".\n\n",
      "Expression cells: ",
      length(cells),
      "\n",
      "Metadata rows: ",
      length(metadata_ids),
      "\n",
      "Expression cells not found in metadata: ",
      length(missing_in_metadata),
      "\n",
      "Metadata barcodes not found in matrix: ",
      length(metadata_not_in_matrix),
      "\n\n",
      "Example expression barcodes not matched:\n",
      paste(
        head(
          missing_in_metadata,
          10
        ),
        collapse = "\n"
      ),
      "\n\n",
      "No automatic/heuristic barcode transformation was performed."
    )
  )
}


# ==============================================================================
# RHAPSODY CANONICAL METADATA
# ==============================================================================

get_vendor_sample_tag <- function(
  seu,
  configured_column = "Sample_Tag"
) {

  md <- seu[[]]

  candidates <- unique(
    c(
      configured_column,
      "Sample_Tag",
      "SampleTag",
      "Sample Tag"
    )
  )

  hit <- intersect(
    candidates,
    colnames(md)
  )

  if (length(hit) == 0) {

    stop(
      paste0(
        "Could not identify Sample Tag column.\n",
        "Metadata columns available: ",
        paste(
          colnames(md),
          collapse = ", "
        )
      )
    )
  }

  as.character(
    md[[hit[1]]]
  )
}


get_vendor_sample_name <- function(
  seu,
  configured_column = "Sample_Name"
) {

  md <- seu[[]]

  candidates <- unique(
    c(
      configured_column,
      "Sample_Name",
      "SampleName",
      "Sample Name"
    )
  )

  hit <- intersect(
    candidates,
    colnames(md)
  )

  if (length(hit) == 0) {
    return(
      rep(
        NA_character_,
        ncol(seu)
      )
    )
  }

  as.character(
    md[[hit[1]]]
  )
}


add_rhapsody_canonical_fields <- function(
  seu,
  sample_tag_column = "Sample_Tag",
  sample_name_column = "Sample_Name"
) {

  vendor_tag <- get_vendor_sample_tag(
    seu,
    configured_column = sample_tag_column
  )

  vendor_name <- get_vendor_sample_name(
    seu,
    configured_column = sample_name_column
  )


  #
  # Preserve vendor values in canonical copies.
  #

  seu$rhapsody_sample_tag <-
    vendor_tag

  seu$rhapsody_sample_name <-
    vendor_name


  #
  # Assignment state.
  #
  # BD may encode Multiplet/Undetermined in Sample_Tag,
  # Sample_Name, or both.
  #

  call_for_status <- vendor_name

  missing_name <- is.na(call_for_status) |
    trimws(call_for_status) == ""

  call_for_status[
    missing_name
  ] <- vendor_tag[
    missing_name
  ]


  status <- rep(
    "assigned",
    length(call_for_status)
  )

  status[
    is.na(call_for_status) |
      trimws(call_for_status) == ""
  ] <- "unknown"

  status[
    grepl(
      "undetermined",
      call_for_status,
      ignore.case = TRUE
    )
  ] <- "undetermined"

  status[
    grepl(
      "multiplet",
      call_for_status,
      ignore.case = TRUE
    )
  ] <- "multiplet"


  seu$rhapsody_assignment_status <-
    status

  seu$rhapsody_is_multiplet <-
    status == "multiplet"

  seu$rhapsody_multiplet_status <-
    ifelse(
      status == "multiplet",
      "multiplet",
      ifelse(
        status == "assigned",
        "not_called_multiplet",
        "not_assessable"
      )
    )

  seu
}


# ==============================================================================
# TARGET METADATA
#
# Core rule:
#
#     physical capture + Sample Tag
#                 ↓
#        exact targets row
#                 ↓
#    biological/project metadata
#
# Multiplet / Undetermined remain present but are NOT assigned biological
# patient/group/treatment metadata.
# ==============================================================================

add_target_metadata <- function(
  seu,
  targets,
  capture_id,
  capture_column,
  tag_column,
  shortname_column
) {

  md <- seu[[]]

  if (
    !"rhapsody_sample_tag" %in%
    colnames(md)
  ) {

    stop(
      "rhapsody_sample_tag is missing before target mapping."
    )
  }


  target_capture <- targets[
    as.character(
      targets[[capture_column]]
    ) == capture_id,
    ,
    drop = FALSE
  ]


  if (nrow(target_capture) == 0) {

    stop(
      "No targets rows found for capture_id='",
      capture_id,
      "' using targets column '",
      capture_column,
      "'."
    )
  }


  if (
    anyDuplicated(
      as.character(
        target_capture[[tag_column]]
      )
    )
  ) {

    stop(
      paste0(
        "Duplicate Sample Tags within capture ",
        capture_id,
        ". Mapping capture + tag must be unique."
      )
    )
  }


  tag_to_row <- setNames(
    seq_len(
      nrow(target_capture)
    ),
    as.character(
      target_capture[[tag_column]]
    )
  )


  cell_tags <- as.character(
    seu$rhapsody_sample_tag
  )

  matched_row <- unname(
    tag_to_row[
      cell_tags
    ]
  )


  assigned_mask <-
    seu$rhapsody_assignment_status == "assigned"

  mapped_mask <-
    assigned_mask &
    !is.na(matched_row)


  # ---------------------------------------------------------------------------
  # Canonical biological sample ID
  # ---------------------------------------------------------------------------

  biological_sample_id <- rep(
    NA_character_,
    ncol(seu)
  )

  biological_sample_id[
    mapped_mask
  ] <- as.character(
    target_capture[
      matched_row[mapped_mask],
      shortname_column,
      drop = TRUE
    ]
  )

  seu$biological_sample_id <-
    biological_sample_id


  # ---------------------------------------------------------------------------
  # Mapping status
  # ---------------------------------------------------------------------------

  mapping_status <- rep(
    "not_applicable",
    ncol(seu)
  )

  mapping_status[
    assigned_mask
  ] <- "unmapped"

  mapping_status[
    mapped_mask
  ] <- "mapped"

  seu$target_mapping_status <-
    mapping_status


  # ---------------------------------------------------------------------------
  # Propagate ALL target columns to mapped singlets.
  #
  # Vendor metadata is never overwritten.
  # If target column already exists, target_<column> is used.
  # ---------------------------------------------------------------------------

  for (
    target_col in colnames(target_capture)
  ) {

    output_col <- target_col

    if (
      output_col %in%
      colnames(seu[[]])
    ) {

      output_col <- paste0(
        "target_",
        target_col
      )
    }


    values <- rep(
      NA,
      ncol(seu)
    )


    source_values <-
      target_capture[[target_col]]


    if (
      is.character(source_values) ||
      is.factor(source_values)
    ) {

      values <- rep(
        NA_character_,
        ncol(seu)
      )

      values[
        mapped_mask
      ] <- as.character(
        source_values[
          matched_row[mapped_mask]
        ]
      )

    } else {

      values <- rep(
        NA_real_,
        ncol(seu)
      )

      values[
        mapped_mask
      ] <- as.numeric(
        source_values[
          matched_row[mapped_mask]
        ]
      )
    }


    seu[[output_col]] <- values
  }


  #
  # Assigned vendor tags that cannot be mapped are preserved,
  # but explicitly flagged.
  #

  unresolved <- assigned_mask &
    is.na(matched_row)

  if (any(unresolved)) {

    warning(
      sum(unresolved),
      " assigned Rhapsody cells could not be mapped to targets for capture ",
      capture_id,
      "."
    )
  }


  seu
}



# ==============================================================================
# UNMAPPED ASSIGNMENT AUDIT
# ==============================================================================

build_unmapped_assignments_summary <- function(seu) {

  md <- seu[[]]

  required_cols <- c(
    "capture_id",
    "rhapsody_sample_tag",
    "rhapsody_sample_name",
    "rhapsody_assignment_status",
    "target_mapping_status"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(md)
  )

  if (length(missing_cols) > 0) {

    stop(
      "Cannot build unmapped assignment audit. Missing metadata columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  unmapped <- md[
    md$rhapsody_assignment_status == "assigned" &
      md$target_mapping_status == "unmapped",
    required_cols,
    drop = FALSE
  ]

  if (nrow(unmapped) == 0) {

    return(
      data.frame(
        capture_id = character(0),
        rhapsody_sample_tag = character(0),
        rhapsody_sample_name = character(0),
        n_cells = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }

  key <- data.frame(
    capture_id = as.character(unmapped$capture_id),
    rhapsody_sample_tag = as.character(unmapped$rhapsody_sample_tag),
    rhapsody_sample_name = as.character(unmapped$rhapsody_sample_name),
    stringsAsFactors = FALSE
  )

  key$rhapsody_sample_name_key <- ifelse(
    is.na(key$rhapsody_sample_name) |
      trimws(key$rhapsody_sample_name) == "",
    "__NA__",
    key$rhapsody_sample_name
  )

  counts <- aggregate(
    rep(1L, nrow(key)),
    by = list(
      capture_id = key$capture_id,
      rhapsody_sample_tag = key$rhapsody_sample_tag,
      rhapsody_sample_name_key = key$rhapsody_sample_name_key
    ),
    FUN = sum
  )

  colnames(counts)[4] <- "n_cells"

  counts$rhapsody_sample_name <- ifelse(
    counts$rhapsody_sample_name_key == "__NA__",
    NA_character_,
    counts$rhapsody_sample_name_key
  )

  counts <- counts[
    ,
    c(
      "capture_id",
      "rhapsody_sample_tag",
      "rhapsody_sample_name",
      "n_cells"
    ),
    drop = FALSE
  ]

  counts <- counts[
    order(
      counts$capture_id,
      counts$rhapsody_sample_tag,
      counts$rhapsody_sample_name,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(counts) <- NULL

  counts
}


# ==============================================================================
# CELL IDENTIFIERS
# ==============================================================================

canonicalize_cell_ids <- function(
  seu,
  capture_id
) {

  original_cells <- colnames(
    seu
  )

  seu$source_barcode <-
    original_cells


  prefix <- paste0(
    capture_id,
    "_"
  )


  if (
    all(
      startsWith(
        original_cells,
        prefix
      )
    )
  ) {

    new_cells <- original_cells

  } else {

    new_cells <- paste0(
      prefix,
      original_cells
    )
  }


  if (
    anyDuplicated(
      new_cells
    )
  ) {

    stop(
      "Duplicated cell IDs after capture prefixing: ",
      capture_id
    )
  }


  seu <- RenameCells(
    object = seu,
    new.names = new_cells
  )


  seu$cell_id <-
    colnames(seu)

  seu
}


# ==============================================================================
# LOAD ONE BD RHAPSODY CAPTURE
# ==============================================================================

load_bd_rhapsody_capture <- function(
  capture_row,
  targets,
  config,
  project_dir
) {

  capture_id <- as.character(
    capture_row$capture_id[[1]]
  )


  message_step(
    paste0(
      "BD Rhapsody capture: ",
      capture_id
    )
  )


  matrix_path <- resolve_path(
    project_dir,
    capture_row$matrix_path[[1]]
  )


  sample_tag_path <- resolve_path(
    project_dir,
    capture_row$sample_tag_calls[[1]]
  )


  barcode_prefix <- NA_character_

  if (
    "barcode_prefix" %in%
    colnames(capture_row)
  ) {

    barcode_prefix <-
      as.character(
        capture_row$barcode_prefix[[1]]
      )
  }


  sample_tag_skip <- 7

  if (
    "sample_tag_skip" %in%
    colnames(capture_row) &&
    !is_missing_value(
      capture_row$sample_tag_skip[[1]]
    )
  ) {

    sample_tag_skip <- as.integer(
      capture_row$sample_tag_skip[[1]]
    )
  }


  cell_type_path <- NA_character_

  if (
    "cell_type_file" %in%
    colnames(capture_row)
  ) {

    cell_type_path <- resolve_path(
      project_dir,
      capture_row$cell_type_file[[1]]
    )
  }


  # ---------------------------------------------------------------------------
  # Primary matrices
  # ---------------------------------------------------------------------------

  primary <- read_primary_matrix(
    matrix_path
  )


  rna_counts <- extract_gene_expression(
    primary
  )


  message(
    "RNA features: ",
    nrow(rna_counts)
  )

  message(
    "Putative cells loaded: ",
    ncol(rna_counts)
  )


  seu <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = capture_id
  )


  # ---------------------------------------------------------------------------
  # Optional AB / ADT
  # ---------------------------------------------------------------------------

  antibody_counts <- extract_antibody_capture(
    primary
  )


  if (!is.null(antibody_counts)) {

    if (
      !identical(
        colnames(antibody_counts),
        colnames(rna_counts)
      )
    ) {

      stop(
        "RNA and Antibody Capture cell barcodes differ for capture ",
        capture_id
      )
    }


    seu[["ADT"]] <- CreateAssayObject(
      counts = antibody_counts
    )


    message(
      "ADT assay loaded: ",
      nrow(antibody_counts),
      " features"
    )
  }


  rm(
    primary,
    rna_counts,
    antibody_counts
  )

  invisible(gc())


  # ---------------------------------------------------------------------------
  # Sample Tag metadata
  # ---------------------------------------------------------------------------

  tags <- read_rhapsody_sample_tags(
    sample_tag_path,
    skip_lines = sample_tag_skip
  )


  tags <- align_metadata_to_cells(
    metadata = tags,
    cells = colnames(seu),
    label = paste0(
      "Sample Tag metadata / ",
      capture_id
    ),
    barcode_prefix = barcode_prefix
  )


  seu <- AddMetaData(
    object = seu,
    metadata = tags
  )


  rm(tags)

  invisible(gc())


  # ---------------------------------------------------------------------------
  # Optional experimental cell-type metadata
  # ---------------------------------------------------------------------------

  if (!is_missing_value(cell_type_path)) {

    cell_type_metadata <-
      read_optional_metadata(
        cell_type_path
      )


    cell_type_metadata <-
      align_metadata_to_cells(
        metadata = cell_type_metadata,
        cells = colnames(seu),
        label = paste0(
          "experimental cell-type metadata / ",
          capture_id
        ),
        barcode_prefix = barcode_prefix
      )


    seu <- AddMetaData(
      object = seu,
      metadata = cell_type_metadata
    )


    rm(
      cell_type_metadata
    )

    invisible(gc())
  }


  # ---------------------------------------------------------------------------
  # Canonical physical/provenance metadata
  # ---------------------------------------------------------------------------

  seu$capture_id <-
    capture_id

  seu$platform <-
    "BD_Rhapsody"

  seu$input_type <-
    "bd_rhapsody"


  # ---------------------------------------------------------------------------
  # Rhapsody canonical evidence
  # ---------------------------------------------------------------------------

  configured_sample_tag_column <-
    config$rhapsody$sample_tag_column %||%
    "Sample_Tag"

  configured_sample_name_column <-
    config$rhapsody$sample_name_column %||%
    "Sample_Name"


  seu <- add_rhapsody_canonical_fields(
    seu,
    sample_tag_column =
      configured_sample_tag_column,
    sample_name_column =
      configured_sample_name_column
  )


  # ---------------------------------------------------------------------------
  # Biological metadata
  # ---------------------------------------------------------------------------

  seu <- add_target_metadata(
    seu = seu,
    targets = targets,
    capture_id = capture_id,
    capture_column =
      config$project$capture_column,
    tag_column =
      config$project$tag_column,
    shortname_column =
      config$project$shortname_column
  )


  # ---------------------------------------------------------------------------
  # Unique identity
  # ---------------------------------------------------------------------------

  seu <- canonicalize_cell_ids(
    seu,
    capture_id
  )


  seu
}


# ==============================================================================
# LEGACY SEURAT RDS SUPPORT
# ==============================================================================

get_counts_matrix <- function(obj) {

  assay_name <- DefaultAssay(
    obj
  )


  tryCatch(
    {

      GetAssayData(
        obj,
        assay = assay_name,
        layer = "counts"
      )

    },
    error = function(e) {

      GetAssayData(
        obj,
        assay = assay_name,
        slot = "counts"
      )
    }
  )
}


find_seurat_rds <- function(
  sample_dir,
  sample_id
) {

  rds_files <- list.files(
    sample_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )


  if (length(rds_files) == 0) {

    stop(
      "No RDS found for ",
      sample_id,
      " in ",
      sample_dir
    )
  }


  if (length(rds_files) > 1) {

    stop(
      paste0(
        "Multiple RDS files found for ",
        sample_id,
        ":\n",
        paste(
          rds_files,
          collapse = "\n"
        )
      )
    )
  }


  rds_files[[1]]
}


make_lightweight_seurat_object <- function(
  raw_obj,
  sample_id
) {

  counts <- get_counts_matrix(
    raw_obj
  )


  original_metadata <-
    raw_obj[[]]


  if (
    !all(
      colnames(counts) %in%
      rownames(original_metadata)
    )
  ) {

    stop(
      "Seurat metadata does not contain all expression barcodes for ",
      sample_id
    )
  }


  original_metadata <-
    original_metadata[
      colnames(counts),
      ,
      drop = FALSE
    ]


  if (
    "orig.ident" %in%
    colnames(original_metadata)
  ) {

    original_metadata$source_orig_ident <-
      original_metadata$orig.ident
  }


  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    meta.data = original_metadata
  )


  obj$input_type <-
    "seurat_rds"

  obj$sample_id <-
    sample_id


  obj <- canonicalize_cell_ids(
    obj,
    sample_id
  )


  obj
}


# ==============================================================================
# MAIN
# ==============================================================================

seu <- NULL

summary_list <- list()

capture_summary_list <- list()


# ==============================================================================
# BD RHAPSODY ROUTE
# ==============================================================================

if (input_type == "bd_rhapsody") {

  capture_manifest_path <- resolve_path(
    project_dir,
    config$input$capture_manifest
  )


  check_file(
    capture_manifest_path,
    "BD Rhapsody capture manifest"
  )


  captures <- readr::read_csv(
    capture_manifest_path,
    show_col_types = FALSE
  )


  required_capture_columns <- c(
    "capture_id",
    "matrix_path",
    "sample_tag_calls"
  )


  missing_capture_columns <- setdiff(
    required_capture_columns,
    colnames(captures)
  )


  if (
    length(
      missing_capture_columns
    ) > 0
  ) {

    stop(
      "Missing capture manifest columns: ",
      paste(
        missing_capture_columns,
        collapse = ", "
      )
    )
  }


  if (
    anyDuplicated(
      captures$capture_id
    )
  ) {

    stop(
      "capture_id must be unique in capture manifest."
    )
  }


  message_step(
    "Loading BD Rhapsody captures"
  )


  message(
    "Number of captures: ",
    nrow(captures)
  )


  for (
    i in seq_len(
      nrow(captures)
    )
  ) {

    capture_row <- captures[
      i,
      ,
      drop = FALSE
    ]


    capture_id <- as.character(
      capture_row$capture_id[[1]]
    )


    obj <- load_bd_rhapsody_capture(
      capture_row = capture_row,
      targets = targets,
      config = config,
      project_dir = project_dir
    )


    n_multiplet <- sum(
      obj$rhapsody_assignment_status ==
        "multiplet",
      na.rm = TRUE
    )


    n_undetermined <- sum(
      obj$rhapsody_assignment_status ==
        "undetermined",
      na.rm = TRUE
    )


    n_mapped <- sum(
      obj$target_mapping_status ==
        "mapped",
      na.rm = TRUE
    )


    n_unmapped <- sum(
      obj$target_mapping_status ==
        "unmapped",
      na.rm = TRUE
    )


    capture_summary_list[[capture_id]] <-
      data.frame(

        capture_id =
          capture_id,

        cells =
          ncol(obj),

        rna_features =
          nrow(obj),

        ab_features =
          if (
            "ADT" %in%
            Assays(obj)
          ) {
            nrow(
              obj[["ADT"]]
            )
          } else {
            0
          },

        mapped_cells =
          n_mapped,

        unmapped_assigned_cells =
          n_unmapped,

        multiplets =
          n_multiplet,

        undetermined =
          n_undetermined,

        stringsAsFactors =
          FALSE
      )


    if (is.null(seu)) {

      seu <- obj

    } else {

      current_assays <- sort(
        Assays(seu)
      )

      incoming_assays <- sort(
        Assays(obj)
      )


      if (
        !identical(
          current_assays,
          incoming_assays
        )
      ) {

        stop(
          paste0(
            "Assay mismatch during capture merge.\n",
            "Current assays: ",
            paste(
              current_assays,
              collapse = ", "
            ),
            "\nIncoming capture ",
            capture_id,
            ": ",
            paste(
              incoming_assays,
              collapse = ", "
            )
          )
        )
      }


      message(
        "Incremental merge: ",
        capture_id
      )


      seu <- merge(
        x = seu,
        y = obj,
        project = project_name,
        merge.data = FALSE,
        merge.dr = FALSE
      )
    }


    rm(obj)

    invisible(gc())


    message(
      "Current merged cells: ",
      ncol(seu)
    )
  }


# ==============================================================================
# LEGACY SEURAT RDS ROUTE
# ==============================================================================

} else if (
  input_type == "seurat_rds"
) {

  count_dir <- resolve_path(
    project_dir,
    config$input$count_dir
  )


  check_dir(
    count_dir,
    "Seurat RDS directory"
  )


  short_column <-
    config$project$shortname_column


  sample_ids <-
    as.character(
      targets[[short_column]]
    )


  if (
    anyDuplicated(
      sample_ids
    )
  ) {

    stop(
      "Duplicate sample identifiers in targets column ",
      short_column
    )
  }


  for (
    sample_id in sample_ids
  ) {

    message_step(
      paste0(
        "Legacy Seurat RDS: ",
        sample_id
      )
    )


    sample_dir <- file.path(
      count_dir,
      sample_id
    )


    check_dir(
      sample_dir,
      paste0(
        "sample directory ",
        sample_id
      )
    )


    rds_path <- find_seurat_rds(
      sample_dir,
      sample_id
    )


    raw_obj <- readRDS(
      rds_path
    )


    obj <- make_lightweight_seurat_object(
      raw_obj,
      sample_id
    )


    rm(raw_obj)

    invisible(gc())


    row_meta <- targets[
      as.character(
        targets[[short_column]]
      ) == sample_id,
      ,
      drop = FALSE
    ]


    if (
      nrow(row_meta) != 1
    ) {

      stop(
        "Expected exactly one targets row for ",
        sample_id
      )
    }


    for (
      column in colnames(row_meta)
    ) {

      output_column <- column


      if (
        output_column %in%
        colnames(obj[[]])
      ) {

        output_column <- paste0(
          "project_",
          column
        )
      }


      value <- row_meta[[column]][1]


      obj[[output_column]] <- rep(
        value,
        ncol(obj)
      )
    }


    summary_list[[sample_id]] <-
      data.frame(

        sample_id =
          sample_id,

        cells =
          ncol(obj),

        genes =
          nrow(obj),

        source =
          rds_path,

        stringsAsFactors =
          FALSE
      )


    if (is.null(seu)) {

      seu <- obj

    } else {

      seu <- merge(
        x = seu,
        y = obj,
        project = project_name,
        merge.data = FALSE,
        merge.dr = FALSE
      )
    }


    rm(obj)

    invisible(gc())
  }


# ==============================================================================
# UNSUPPORTED INPUT
# ==============================================================================

} else {

  stop(
    "Unsupported input.type: ",
    input_type,
    ". Supported values: bd_rhapsody, seurat_rds"
  )
}


# ==============================================================================
# FINAL VALIDATION
# ==============================================================================

if (is.null(seu)) {

  stop(
    "No cells were loaded."
  )
}


if (
  ncol(seu) == 0
) {

  stop(
    "Final Seurat object contains zero cells."
  )
}


if (
  anyDuplicated(
    colnames(seu)
  )
) {

  stop(
    "Duplicated final cell identifiers detected."
  )
}


if (
  !"cell_id" %in%
  colnames(seu[[]])
) {

  stop(
    "Canonical cell_id metadata is missing."
  )
}


if (
  anyDuplicated(
    seu$cell_id
  )
) {

  stop(
    "Duplicated cell_id values detected."
  )
}


message_step(
  "Final canonical object"
)


message(
  "Cells: ",
  ncol(seu)
)


message(
  "RNA features: ",
  nrow(seu)
)


message(
  "Assays: ",
  paste(
    Assays(seu),
    collapse = ", "
  )
)


message(
  "Metadata columns: ",
  ncol(
    seu[[]]
  )
)


# ==============================================================================
# RHAPSODY FINAL AUDIT
# ==============================================================================

if (
  input_type == "bd_rhapsody"
) {

  required_rhapsody_metadata <- c(
    "capture_id",
    "rhapsody_sample_tag",
    "rhapsody_assignment_status",
    "rhapsody_is_multiplet",
    "biological_sample_id",
    "target_mapping_status",
    "source_barcode",
    "cell_id"
  )


  missing_rhapsody_metadata <- setdiff(
    required_rhapsody_metadata,
    colnames(
      seu[[]]
    )
  )


  if (
    length(
      missing_rhapsody_metadata
    ) > 0
  ) {

    stop(
      "Missing required canonical Rhapsody metadata: ",
      paste(
        missing_rhapsody_metadata,
        collapse = ", "
      )
    )
  }


  message(
    "Mapped biological cells: ",
    sum(
      seu$target_mapping_status ==
        "mapped",
      na.rm = TRUE
    )
  )


  message(
    "Rhapsody Multiplets preserved: ",
    sum(
      seu$rhapsody_assignment_status ==
        "multiplet",
      na.rm = TRUE
    )
  )


  message(
    "Rhapsody Undetermined preserved: ",
    sum(
      seu$rhapsody_assignment_status ==
        "undetermined",
      na.rm = TRUE
    )
  )
}


# ==============================================================================
# SAVE RDS
# ==============================================================================

message_step(
  "Saving canonical Seurat object"
)


saveRDS(
  seu,
  output_rds
)


message(
  "Saved: ",
  output_rds
)


# ==============================================================================
# MAIN SUMMARY
# ==============================================================================

if (
  input_type == "bd_rhapsody"
) {

  main_summary <- data.frame(

    input_type =
      input_type,

    project =
      project_name,

    cells =
      ncol(seu),

    rna_features =
      nrow(seu),

    captures =
      length(
        unique(
          seu$capture_id
        )
      ),

    multiplets =
      sum(
        seu$rhapsody_assignment_status ==
          "multiplet",
        na.rm = TRUE
      ),

    undetermined =
      sum(
        seu$rhapsody_assignment_status ==
          "undetermined",
        na.rm = TRUE
      ),

    mapped_cells =
      sum(
        seu$target_mapping_status ==
          "mapped",
        na.rm = TRUE
      ),

    unmapped_assigned_cells =
      sum(
        seu$target_mapping_status ==
          "unmapped",
        na.rm = TRUE
      ),

    stringsAsFactors =
      FALSE
  )


  readr::write_csv(
    main_summary,
    summary_file
  )


} else {

  readr::write_csv(
    dplyr::bind_rows(
      summary_list
    ),
    summary_file
  )
}


# ==============================================================================
# CAPTURE SUMMARY
# ==============================================================================

if (
  length(
    capture_summary_list
  ) > 0
) {

  capture_summary <- dplyr::bind_rows(
    capture_summary_list
  )


  readr::write_csv(
    capture_summary,
    capture_summary_file
  )
}


# ==============================================================================
# RHAPSODY ASSIGNMENT SUMMARY
# ==============================================================================

if (
  input_type == "bd_rhapsody"
) {

  assignment_summary <- as.data.frame(
    table(
      capture_id =
        seu$capture_id,
      assignment =
        seu$rhapsody_assignment_status,
      useNA =
        "ifany"
    )
  )


  readr::write_csv(
    assignment_summary,
    assignment_summary_file
  )
}



# ==============================================================================
# UNMAPPED ASSIGNMENT SUMMARY
# ==============================================================================

if (
  input_type == "bd_rhapsody"
) {

  unmapped_assignments_summary <-
    build_unmapped_assignments_summary(
      seu
    )


  readr::write_csv(
    unmapped_assignments_summary,
    unmapped_assignments_file,
    na = ""
  )


  message(
    "Unmapped assigned cells summary written: ",
    unmapped_assignments_file,
    " (",
    sum(
      unmapped_assignments_summary$n_cells
    ),
    " cells)"
  )
}


# ==============================================================================
# SESSION INFO
# ==============================================================================

sink(
  session_file
)


cat(
  "project=",
  project_name,
  "\n",
  sep = ""
)


cat(
  "input_type=",
  input_type,
  "\n",
  sep = ""
)


cat(
  "cells=",
  ncol(seu),
  "\n",
  sep = ""
)


cat(
  "features=",
  nrow(seu),
  "\n",
  sep = ""
)


cat(
  "\n"
)


print(
  sessionInfo()
)


sink()


message_step(
  "00_load_data finished successfully"
)
