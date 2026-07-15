#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ============================================================
# 01_build_drs_inputs.R
#
# Purpose:
#   Build the 3 required DRS input contracts from the current IMRS project:
#     1) 01_dataset_manifest.tsv
#     2) 02_sample_metadata.tsv
#     3) 03_expression_matrix_index.tsv
#   and copy/link per-dataset expression matrices into:
#     DRS/inputs/expression_matrices/<dataset_id>.tsv
#
# Prototype assumptions:
#   - source counts are raw counts from Step 3
#   - one scoring design file per dataset:
#       00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#   - one count matrix per dataset:
#       03_counts/<DATASET>/**/gene_counts_clean.tsv
#
# Notes:
#   - this script is intentionally permissive and heuristic
#   - it tries to infer group / replicate / time / tissue / organism
#   - if a field is missing, it fills with "unknown" or NA
#
# DRS requires these contracts and raw-count declarations for portability,
# with per-dataset expression matrices referenced by an index. :contentReference[oaicite:0]{index=0}
# The IMRS Step 3 output is exactly the gene-level count matrix needed here. 
# ============================================================

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

design_root <- file.path(project_root, "00_metadata", "verified_metadata", "scoring")
counts_root <- file.path(project_root, "03_counts")

drs_root  <- file.path(project_root, "DRS")
input_dir <- file.path(drs_root, "inputs")
expr_dir  <- file.path(input_dir, "expression_matrices")
log_dir   <- file.path(drs_root, "logs")

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

manifest_out <- file.path(input_dir, "01_dataset_manifest.tsv")
sample_out   <- file.path(input_dir, "02_sample_metadata.tsv")
index_out    <- file.path(input_dir, "03_expression_matrix_index.tsv")
issues_out   <- file.path(log_dir, "01_build_drs_inputs_issues.tsv")

# Set TRUE to copy count files into DRS/inputs/expression_matrices.
# Set FALSE to just reference original count files directly in the index.
copy_expression_matrices <- TRUE

# -------------------------
# HELPERS
# -------------------------
read_table_robust <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)

  x <- tryCatch(read_tsv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 2) return(x)

  x <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 2) return(x)

  x <- tryCatch(read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 2) return(x)

  stop("Could not read table: ", path)
}

safe_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

as_bool_string <- function(x) {
  ifelse(isTRUE(x), "TRUE", "FALSE")
}

norm_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

infer_dataset_class <- function(dataset_id) {
  case_when(
    dataset_id %in% c("GSE39129", "GSE167521", "GSE264344") ~ "anchor",
    dataset_id %in% c("GSE279744") ~ "calibration",
    dataset_id %in% c("GSE190850") ~ "transfer",
    TRUE ~ "external"
  )
}

infer_assay_type <- function(dataset_id, design_df = NULL) {
  "bulk_rnaseq"
}

infer_group_column <- function(design_df) {
  safe_col(design_df, c(
    "group",
    "condition_simple",
    "condition",
    "treatment_group",
    "group_label"
  ))
}

infer_group_values <- function(design_df) {
  grp_col <- infer_group_column(design_df)
  if (is.null(grp_col)) return(rep(NA_character_, nrow(design_df)))

  x <- toupper(norm_text(design_df[[grp_col]]))

  x <- case_when(
    x %in% c("CONTROL", "CTRL", "UNTREATED", "NAIVE", "PBS", "VEHICLE", "BASELINE", "PRE") ~ "CONTROL",
    x %in% c("DELIVERY", "TREATED", "TREATMENT", "CASE", "POST") ~ "DELIVERY",
    TRUE ~ x
  )
  x
}

infer_replicate <- function(design_df, group_vec) {
  rep_col <- safe_col(design_df, c("replicate", "rep", "biological_replicate"))
  if (!is.null(rep_col)) {
    out <- suppressWarnings(as.integer(design_df[[rep_col]]))
    if (!all(is.na(out))) return(out)
  }

  # fallback: rank within group
  tibble(group = group_vec) %>%
    group_by(group) %>%
    mutate(replicate = row_number()) %>%
    ungroup() %>%
    pull(replicate)
}

infer_time_hr <- function(design_df) {
  time_col <- safe_col(design_df, c(
    "time_hr", "timepoint_hr", "hours", "hour", "timepoint", "time_h", "time"
  ))
  if (is.null(time_col)) return(rep(NA_real_, nrow(design_df)))

  x <- design_df[[time_col]]

  if (is.numeric(x)) return(as.numeric(x))

  x <- as.character(x)
  x <- str_replace_all(x, "(?i)hours|hour|hrs|hr|h", "")
  x <- str_trim(x)
  suppressWarnings(as.numeric(x))
}

infer_batch <- function(design_df) {
  batch_col <- safe_col(design_df, c("batch", "batch_id", "lane", "run", "seq_batch"))
  if (is.null(batch_col)) return(rep(NA_character_, nrow(design_df)))
  norm_text(design_df[[batch_col]])
}

infer_sample_id <- function(design_df) {
  sid_col <- safe_col(design_df, c("sample_id", "sample", "run_id", "srr", "SRR"))
  if (is.null(sid_col)) stop("Could not find sample_id-like column in design file.")
  norm_text(design_df[[sid_col]])
}

infer_dataset_level_field <- function(design_df, candidates, default = "unknown") {
  col <- safe_col(design_df, candidates)
  if (is.null(col)) return(default)

  vals <- norm_text(design_df[[col]])
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(default)

  # use first unique non-missing value if consistent, else "mixed"
  u <- unique(vals)
  if (length(u) == 1) return(u[1])
  "mixed"
}

infer_organism <- function(dataset_id, design_df) {
  x <- infer_dataset_level_field(design_df, c("organism", "species"), default = NA_character_)

  if (!is.na(x)) {
    x2 <- tolower(x)
    if (str_detect(x2, "human|homo")) return("human")
    if (str_detect(x2, "mouse|mus")) return("mouse")
    if (str_detect(x2, "rat")) return("rat")
    return(x2)
  }

  if (dataset_id %in% c("GSE190850")) return("human")
  "mouse"
}

infer_tissue <- function(design_df) {
  infer_dataset_level_field(design_df, c("tissue", "organ", "source_name", "sample_type"), default = "unknown")
}

infer_treatment_type <- function(design_df, dataset_id) {
  x <- infer_dataset_level_field(design_df, c("delivery_type", "treatment_type", "treatment", "vector", "platform"),
                                 default = NA_character_)
  if (!is.na(x)) return(x)

  case_when(
    dataset_id == "GSE39129"  ~ "lentiviral_vector",
    dataset_id == "GSE167521" ~ "lnp",
    dataset_id == "GSE264344" ~ "adenovirus",
    dataset_id == "GSE279744" ~ "lnp_mrna",
    dataset_id == "GSE190850" ~ "adenovirus_vaccine",
    TRUE ~ "unknown"
  )
}

find_best_counts_file <- function(dataset_id) {
  root <- file.path(counts_root, dataset_id)
  if (!dir.exists(root)) return(NA_character_)

  files <- list.files(root, pattern = "^gene_counts_clean\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)

  # Prefer validation folder if present, otherwise shortest path.
  score <- ifelse(str_detect(files, "/validation/|\\\\validation\\\\"), 0, 1)
  files[order(score, nchar(files))][1]
}

read_count_dims <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE, n_max = 5)
  full <- read_tsv(path, show_col_types = FALSE, progress = FALSE)

  if (ncol(full) < 2) stop("Count matrix has <2 columns: ", path)

  list(
    n_genes = nrow(full),
    n_samples = ncol(full) - 1L,
    sample_cols = colnames(full)[-1]
  )
}

make_relative_to_project <- function(path_abs) {
  rel <- sub(
    paste0("^", gsub("\\\\", "/", normalizePath(project_root, winslash = "/", mustWork = TRUE)), "/?"),
    "",
    gsub("\\\\", "/", normalizePath(path_abs, winslash = "/", mustWork = TRUE))
  )
  rel
}

# -------------------------
# DISCOVER DATASETS
# -------------------------
design_dirs <- list.dirs(design_root, full.names = TRUE, recursive = FALSE)
dataset_ids <- basename(design_dirs)

if (length(dataset_ids) == 0) {
  stop("No dataset folders found under: ", design_root)
}

message("Found datasets: ", paste(dataset_ids, collapse = ", "))

manifest_rows <- list()
sample_rows   <- list()
index_rows    <- list()
issue_rows    <- list()

# -------------------------
# MAIN LOOP
# -------------------------
for (dataset_id in dataset_ids) {
  message("Processing ", dataset_id, " ...")

  design_path <- file.path(design_root, dataset_id, paste0(dataset_id, "_design.tsv"))
  counts_path <- find_best_counts_file(dataset_id)

  if (!file.exists(design_path)) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "missing_design_file",
      detail = design_path
    )
    next
  }

  if (is.na(counts_path) || !file.exists(counts_path)) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "missing_counts_file",
      detail = as.character(counts_path)
    )
    next
  }

  design_df <- tryCatch(read_table_robust(design_path), error = function(e) e)
  if (inherits(design_df, "error")) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "design_read_failed",
      detail = conditionMessage(design_df)
    )
    next
  }

  # infer fields
  sample_id  <- tryCatch(infer_sample_id(design_df), error = function(e) NULL)
  group_vec  <- infer_group_values(design_df)
  replicate  <- infer_replicate(design_df, group_vec)
  batch_vec  <- infer_batch(design_df)
  time_hr    <- infer_time_hr(design_df)

  if (is.null(sample_id)) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "missing_sample_id_column",
      detail = design_path
    )
    next
  }

  organism       <- infer_organism(dataset_id, design_df)
  assay_type     <- infer_assay_type(dataset_id, design_df)
  tissue         <- infer_tissue(design_df)
  treatment_type <- infer_treatment_type(design_df, dataset_id)
  dataset_class  <- infer_dataset_class(dataset_id)

  # counts info
  count_info <- tryCatch(read_count_dims(counts_path), error = function(e) e)
  if (inherits(count_info, "error")) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "counts_read_failed",
      detail = conditionMessage(count_info)
    )
    next
  }

  count_samples <- norm_text(count_info$sample_cols)

  # Keep only design samples that appear in count matrix
  keep <- sample_id %in% count_samples
  if (!all(keep)) {
    missing_samples <- paste(sample_id[!keep], collapse = ";")
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "design_samples_missing_from_counts",
      detail = missing_samples
    )
  }

  if (sum(keep) == 0) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "no_overlapping_samples",
      detail = counts_path
    )
    next
  }

  sample_id <- sample_id[keep]
  group_vec <- group_vec[keep]
  replicate <- replicate[keep]
  batch_vec <- batch_vec[keep]
  time_hr   <- time_hr[keep]

  # Output expression matrix path
  expr_out_abs <- file.path(expr_dir, paste0(dataset_id, ".tsv"))

  if (copy_expression_matrices) {
    ok <- file.copy(counts_path, expr_out_abs, overwrite = TRUE)
    if (!ok) {
      issue_rows[[length(issue_rows) + 1]] <- tibble(
        dataset_id = dataset_id,
        level = "dataset",
        issue = "failed_to_copy_expression_matrix",
        detail = paste(counts_path, "->", expr_out_abs)
      )
      next
    }
    expr_rel <- make_relative_to_project(expr_out_abs)
  } else {
    expr_rel <- make_relative_to_project(counts_path)
  }

  # manifest row
  manifest_rows[[length(manifest_rows) + 1]] <- tibble(
    dataset_id = dataset_id,
    dataset_class = dataset_class,
    organism = organism,
    assay_type = assay_type,
    tissue = tissue,
    treatment_type = treatment_type,
    grouping_variable = "group",
    include_in_analysis = "TRUE",
    notes = NA_character_
  )

  # sample metadata rows
  sample_rows[[length(sample_rows) + 1]] <- tibble(
    dataset_id = dataset_id,
    sample_id = sample_id,
    group = group_vec,
    replicate = as.integer(replicate),
    batch = batch_vec,
    time_hr = as.numeric(time_hr),
    covariates = NA_character_,
    exclude_sample = "FALSE",
    notes = NA_character_
  )

  # expression index row
  index_rows[[length(index_rows) + 1]] <- tibble(
    dataset_id = dataset_id,
    expression_matrix_file = expr_rel,
    value_type = "raw_counts",
    normalization_method = "none",
    gene_id_type = "as_is",
    n_genes = as.integer(count_info$n_genes),
    n_samples = as.integer(sum(keep)),
    notes = paste0("source_counts_file=", make_relative_to_project(counts_path))
  )

  # basic issues / warnings
  if (any(is.na(group_vec))) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "sample",
      issue = "missing_group_values",
      detail = paste(sample_id[is.na(group_vec)], collapse = ";")
    )
  }

  grp_tab <- table(group_vec, useNA = "no")
  if (any(grp_tab < 2)) {
    issue_rows[[length(issue_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      level = "dataset",
      issue = "small_group_size",
      detail = paste(names(grp_tab), grp_tab, sep = ":", collapse = ";")
    )
  }
}

# -------------------------
# BIND + CLEAN
# -------------------------
manifest_tbl <- bind_rows(manifest_rows) %>%
  distinct(dataset_id, .keep_all = TRUE) %>%
  arrange(dataset_id)

sample_tbl <- bind_rows(sample_rows) %>%
  arrange(dataset_id, group, replicate, sample_id)

index_tbl <- bind_rows(index_rows) %>%
  distinct(dataset_id, .keep_all = TRUE) %>%
  arrange(dataset_id)

issues_tbl <- bind_rows(issue_rows) %>%
  arrange(dataset_id, level, issue)

# -------------------------
# WRITE
# -------------------------
write_tsv(manifest_tbl, manifest_out, na = "")
write_tsv(sample_tbl, sample_out, na = "")
write_tsv(index_tbl, index_out, na = "")
write_tsv(issues_tbl, issues_out, na = "")

message("Saved:")
message("  ", manifest_out)
message("  ", sample_out)
message("  ", index_out)
message("  ", issues_out)

message("")
message("Counts:")
message("  datasets in manifest: ", nrow(manifest_tbl))
message("  rows in sample metadata: ", nrow(sample_tbl))
message("  rows in matrix index: ", nrow(index_tbl))
message("Done.")