#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 08H — Score HUMAN samples using transferred
# mouse->human ortholog weights
#
# built on: scoring.pdf (2026-03-03)
#
# Purpose:
#   External cross-species transfer test only.
#   NO retraining, NO recalibration, NO weight updates.
#
# Inputs:
#   1) Raw integer counts:
#        <project_root>/03_counts/<DATASET>/.../gene_counts_clean.tsv
#   2) Design table with:
#        sample_id
#        condition_simple in {CONTROL, DELIVERY}
#        <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#   3) Human transferred weights:
#        <project_root>/05_score/human_transfer/human_gene_weights.tsv
#   4) Ortholog coverage table:
#        <project_root>/05_score/human_transfer/ortholog_coverage.tsv
#
# Fixed logic:
#   - Normalize counts with DESeq2 size factors on ALL samples, design ~ 1
#   - Transform: log2(norm + 1)
#   - Gene z-score using CONTROLS only
#   - Minimum controls: n >= 3 (else fail)
#   - Gene SD floor: 0.10
#   - Final score SD floor: 0.10
#   - Strip Ensembl version suffix: TRUE
#
# IMPORTANT HUMAN TRANSFER RULE:
#   - Do NOT use mouse 0.80 hard-fail coverage rule here.
#   - Report ortholog coverage and actual scoring overlap.
#   - Interpret conservatively.
#
# Outputs (under 05_score/human_transfer/):
#   - scores/<ID>__imrs_scores.tsv
#   - qc/<ID>__qc_summary.tsv
#   - qc/<ID>__top_contributors.tsv
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# Pass exact dataset folder name as second arg if desired, e.g. GSE190850_HUMAN
dataset_filter <- if (length(args) >= 2) args[2] else NA_character_

weights_path <- file.path(project_root, "05_score", "human_transfer", "human_gene_weights.tsv")
ortholog_cov_path <- file.path(project_root, "05_score", "human_transfer", "ortholog_coverage.tsv")
design_root <- file.path(project_root, "00_metadata", "verified_metadata", "scoring")

min_ctrl_n       <- 3
gene_sd_floor    <- 0.10
score_sd_floor   <- 0.10
strip_ens_ver    <- TRUE
top_k_contrib    <- 20

# Human transfer output base
out_base <- file.path(project_root, "05_score", "human_transfer")
out_scores_dir <- file.path(out_base, "scores")
out_qc_dir     <- file.path(out_base, "qc")
dir.create(out_scores_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_qc_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Helpers
# -------------------------
strip_ensembl_version <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+$", "", x)
}

read_table_robust <- function(path) {
  df <- tryCatch(
    read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(
    read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(
    read_table(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && ncol(df) >= 2) return(df)

  stop("Failed to read table: ", path)
}

find_counts_files <- function(dataset_id) {
  root <- file.path(project_root, "03_counts", dataset_id)
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = "^gene_counts_clean\\.tsv$", recursive = TRUE, full.names = TRUE)
}

make_counts_id <- function(dataset_id, counts_path) {
  root <- file.path(project_root, "03_counts", dataset_id)
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", gsub("\\\\", "/", counts_path))
  rel <- sub("/gene_counts_clean\\.tsv$", "", rel)
  rel <- ifelse(nchar(rel) == 0, "counts", rel)
  rel <- gsub("[^A-Za-z0-9_\\-]+", "_", rel)
  paste0(dataset_id, "__", rel)
}

read_counts_matrix <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (ncol(df) < 2) stop("Counts file has <2 columns: ", path)

  names(df)[1] <- "gene_id"
  gene_id <- as.character(df$gene_id)
  mat <- as.matrix(df[, -1, drop = FALSE])

  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- gene_id
  colnames(mat) <- trimws(colnames(mat))
  mat
}

is_integerish_matrix <- function(mat) {
  all(is.finite(mat)) && all(abs(mat - round(mat)) < 1e-6)
}

load_human_transfer_weights <- function(path) {
  if (!file.exists(path)) stop("Missing human transferred weights: ", path)

  w <- read_tsv(path, show_col_types = FALSE)

  gene_col <- NULL
  for (cand in c("human_gene_id", "gene_id", "gene")) {
    if (cand %in% names(w)) {
      gene_col <- cand
      break
    }
  }
  if (is.null(gene_col)) {
    stop("human_gene_weights.tsv must contain one of: human_gene_id, gene_id, gene")
  }

  weight_col <- NULL
  for (cand in c("weight", "beta_meta")) {
    if (cand %in% names(w)) {
      weight_col <- cand
      break
    }
  }
  if (is.null(weight_col)) {
    stop("human_gene_weights.tsv must contain one of: weight, beta_meta")
  }

  out <- w %>%
    transmute(
      gene_id = as.character(.data[[gene_col]]),
      beta_meta = suppressWarnings(as.numeric(.data[[weight_col]]))
    ) %>%
    filter(!is.na(beta_meta))

  if (strip_ens_ver) {
    out <- out %>% mutate(gene_id = strip_ensembl_version(gene_id))
  }

  out <- out %>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(beta_meta), n = 1, with_ties = FALSE) %>%
    ungroup()

  if (nrow(out) == 0) stop("No usable human transferred weights after cleaning: ", path)
  out
}

load_ortholog_coverage <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  x <- read_tsv(path, show_col_types = FALSE)
  if (!("ortholog_coverage" %in% names(x))) return(NA_real_)
  as.numeric(x$ortholog_coverage[1])
}

# Auto-detect likely human dataset folder if not explicitly provided
detect_human_dataset <- function(counts_root) {
  ds <- list.dirs(counts_root, recursive = FALSE, full.names = FALSE)
  ds <- ds[grepl("^GSE", ds)]

  # Prefer names containing HUMAN
  human_named <- ds[grepl("HUMAN", ds, ignore.case = TRUE)]
  if (length(human_named) == 1) return(human_named)

  # Fallback: prefer exact known transfer dataset folder names if present
  preferred <- c("GSE190850_HUMAN", "GSE190850")
  hit <- preferred[preferred %in% ds]
  if (length(hit) >= 1) return(hit[1])

  NA_character_
}

# -------------------------
# Core scoring function
# -------------------------
score_one_counts_file <- function(dataset_id, counts_path, design_df, weights_df, ortholog_coverage_fixed) {
  counts_id <- make_counts_id(dataset_id, counts_path)

  counts <- read_counts_matrix(counts_path)

  if (strip_ens_ver) {
    rownames(counts) <- strip_ensembl_version(rownames(counts))
  }

  if (!("sample_id" %in% names(design_df))) {
    stop("Design missing sample_id for dataset: ", dataset_id)
  }
  if (!("condition_simple" %in% names(design_df))) {
    stop("Design missing condition_simple (CONTROL/DELIVERY) for dataset: ", dataset_id)
  }

  design <- design_df %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    ) %>%
    filter(sample_id %in% colnames(counts)) %>%
    filter(condition_simple %in% c("CONTROL", "DELIVERY"))

  if (nrow(design) == 0) {
    message("Skip (no overlapping samples): ", counts_id)
    return(invisible(NULL))
  }

  counts <- counts[, design$sample_id, drop = FALSE]

  if (!is_integerish_matrix(counts)) {
    message("FAIL (non-integer counts): ", counts_id)
    qc <- tibble(
      id = counts_id,
      dataset_id = dataset_id,
      counts_path = counts_path,
      pass = FALSE,
      fail_reason = "Counts not integer"
    )
    write_tsv(qc, file.path(out_qc_dir, paste0(counts_id, "__qc_summary.tsv")))
    return(invisible(NULL))
  }

  counts_int <- round(counts)
  storage.mode(counts_int) <- "integer"

  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  if (length(ctrl_ids) < min_ctrl_n) {
    message("FAIL (controls < ", min_ctrl_n, "): ", counts_id, " | n_ctrl=", length(ctrl_ids))
    qc <- tibble(
      id = counts_id,
      dataset_id = dataset_id,
      counts_path = counts_path,
      pass = FALSE,
      fail_reason = paste0("Controls < ", min_ctrl_n),
      n_controls = length(ctrl_ids)
    )
    write_tsv(qc, file.path(out_qc_dir, paste0(counts_id, "__qc_summary.tsv")))
    return(invisible(NULL))
  }

  # Same normalization logic as mouse Step 08
  sf <- estimateSizeFactorsForMatrix(counts_int)
  norm <- sweep(counts_int, 2, sf, "/")
  x <- log2(norm + 1)

  w_genes <- weights_df$gene_id
  overlap <- intersect(rownames(x), w_genes)
  scoring_coverage <- length(overlap) / length(w_genes)

  if (length(overlap) == 0) {
    message("FAIL (no overlapping genes): ", counts_id)
    qc <- tibble(
      id = counts_id,
      dataset_id = dataset_id,
      counts_path = counts_path,
      pass = FALSE,
      fail_reason = "No overlapping genes between human weights and counts",
      ortholog_coverage_fixed = ortholog_coverage_fixed,
      n_transferred_weights = length(w_genes),
      n_overlap = 0,
      scoring_coverage = 0,
      n_controls = length(ctrl_ids)
    )
    write_tsv(qc, file.path(out_qc_dir, paste0(counts_id, "__qc_summary.tsv")))
    return(invisible(NULL))
  }

  wsub <- weights_df %>%
    filter(gene_id %in% overlap) %>%
    arrange(match(gene_id, overlap))

  xsub <- x[overlap, , drop = FALSE]

  ctrl_mat <- xsub[, ctrl_ids, drop = FALSE]
  mu <- rowMeans(ctrl_mat, na.rm = TRUE)
  sdv <- apply(ctrl_mat, 1, sd, na.rm = TRUE)

  sdv_floor <- pmax(sdv, gene_sd_floor)
  n_hit_floor <- sum(sdv < gene_sd_floor, na.rm = TRUE)

  z <- sweep(xsub, 1, mu, "-")
  z <- sweep(z, 1, sdv_floor, "/")

  wvec <- wsub$beta_meta
  names(wvec) <- wsub$gene_id
  wvec <- wvec[rownames(z)]

  imrs_raw <- as.numeric(crossprod(wvec, z))
  names(imrs_raw) <- colnames(z)

  ctrl_scores <- imrs_raw[ctrl_ids]
  mu_score <- mean(ctrl_scores, na.rm = TRUE)
  sd_score <- sd(ctrl_scores, na.rm = TRUE)
  sd_score_floor <- max(sd_score, score_sd_floor)

  imrs_z <- (imrs_raw - mu_score) / sd_score_floor

  out_scores <- tibble(
    sample_id = design$sample_id,
    condition_simple = design$condition_simple,
    imrs_raw = as.numeric(imrs_raw[design$sample_id]),
    imrs_z = as.numeric(imrs_z[design$sample_id]),
    ortholog_coverage_fixed = ortholog_coverage_fixed,
    scoring_coverage = scoring_coverage
  )

  write_tsv(out_scores, file.path(out_scores_dir, paste0(counts_id, "__imrs_scores.tsv")))

  qc <- tibble(
    id = counts_id,
    dataset_id = dataset_id,
    counts_path = counts_path,
    pass = TRUE,
    n_samples = nrow(design),
    n_controls = length(ctrl_ids),
    n_delivery = sum(design$condition_simple == "DELIVERY"),
    n_transferred_weights = length(w_genes),
    n_overlap = length(overlap),
    ortholog_coverage_fixed = ortholog_coverage_fixed,
    scoring_coverage = scoring_coverage,
    n_genes_sd_floor = n_hit_floor,
    gene_sd_floor = gene_sd_floor,
    score_sd_floor = score_sd_floor,
    controls_imrs_z_mean = mean(imrs_z[ctrl_ids], na.rm = TRUE),
    controls_imrs_z_sd = sd(imrs_z[ctrl_ids], na.rm = TRUE)
  )

  write_tsv(qc, file.path(out_qc_dir, paste0(counts_id, "__qc_summary.tsv")))

  contrib_mat <- sweep(z, 1, wvec, "*")
  top_tbl <- map_dfr(colnames(contrib_mat), function(sid) {
    v <- contrib_mat[, sid]
    ord <- order(abs(v), decreasing = TRUE)
    ord <- ord[seq_len(min(top_k_contrib, length(ord)))]
    tibble(
      sample_id = sid,
      gene_id = rownames(contrib_mat)[ord],
      weight = as.numeric(wvec[rownames(contrib_mat)[ord]]),
      z = as.numeric(z[rownames(contrib_mat)[ord], sid]),
      w_times_z = as.numeric(v[ord]),
      abs_w_times_z = abs(as.numeric(v[ord])),
      rank = seq_along(ord)
    )
  })

  write_tsv(top_tbl, file.path(out_qc_dir, paste0(counts_id, "__top_contributors.tsv")))

  message(
    "OK: ", counts_id,
    " | n_samples=", nrow(design),
    " | n_ctrl=", length(ctrl_ids),
    " | ortholog_coverage=", ifelse(is.na(ortholog_coverage_fixed), "NA", sprintf("%.3f", ortholog_coverage_fixed)),
    " | scoring_coverage=", sprintf("%.3f", scoring_coverage),
    " | sd_floor_genes=", n_hit_floor
  )

  invisible(NULL)
}

# -------------------------
# MAIN
# -------------------------
if (!file.exists(weights_path)) stop("Human transferred weights missing: ", weights_path)

weights_df <- load_human_transfer_weights(weights_path)
ortholog_coverage_fixed <- load_ortholog_coverage(ortholog_cov_path)

message(
  "Loaded human transferred weights: n_genes=", nrow(weights_df),
  " | path=", weights_path
)
message(
  "Loaded ortholog coverage: ",
  ifelse(is.na(ortholog_coverage_fixed), "NA", sprintf("%.4f", ortholog_coverage_fixed))
)

counts_root <- file.path(project_root, "03_counts")
if (!dir.exists(counts_root)) stop("Missing 03_counts root: ", counts_root)

if (!dir.exists(design_root)) {
  stop("Missing design_root: ", design_root)
}

# Restrict to one human dataset only
if (is.na(dataset_filter)) {
  dataset_filter <- detect_human_dataset(counts_root)
}

if (is.na(dataset_filter) || !nzchar(dataset_filter)) {
  stop("Could not determine human dataset automatically. Please pass it as second argument, e.g. GSE190850_HUMAN")
}

dataset_ids <- dataset_filter

message("Human dataset selected: ", dataset_ids)

for (ds in dataset_ids) {
  design_path <- file.path(design_root, ds, paste0(ds, "_design.tsv"))

  if (!file.exists(design_path)) {
    stop("Missing scoring design for selected dataset: ", ds, " | ", design_path)
  }

  design_df <- read_table_robust(design_path)

  if (!("condition_simple" %in% names(design_df))) {
    stop("Design missing condition_simple for selected dataset: ", ds)
  }

  cfiles <- find_counts_files(ds)
  if (length(cfiles) == 0) {
    stop("No gene_counts_clean.tsv found for selected dataset: ", ds)
  }

  message("\n==============================")
  message("DATASET: ", ds, " | count_files=", length(cfiles))
  message("==============================")

  for (cf in cfiles) {
    tryCatch(
      score_one_counts_file(ds, cf, design_df, weights_df, ortholog_coverage_fixed),
      error = function(e) {
        message("ERROR scoring: ", ds, " | ", cf, " | ", conditionMessage(e))
      }
    )
  }
}

message("\nDONE: Step 08H human transfer scoring finished (strict transfer; no retraining).")