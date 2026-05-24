#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 08 — Score samples using FROZEN ANCHOR WEIGHTS
#
# built on: scoring.pdf (2026-03-03)
#
# Inputs (per dataset):
#   1) Raw integer counts:
#        <project_root>/03_counts/<DATASET>/.../gene_counts_clean.tsv
#   2) Design table with:
#        sample_id
#        condition_simple  in {CONTROL, DELIVERY}
#      (UPDATED PATH)
#        <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#   3) Frozen anchor weights (NEVER refit):
#        <project_root>/05_score/anchors/gene_weights.tsv
#
# Fixed defaults:
#   - Normalize counts with DESeq2 size factors on ALL samples, design ~ 1
#   - Transform: log2(norm + 1)
#   - Gene z-score using CONTROLS only
#   - Minimum controls: n >= 3 (else fail)
#   - Gene SD floor: 0.10
#   - Final score SD floor: 0.10
#   - Strip Ensembl version suffix: TRUE
#   - Coverage >= 0.80 (else fail)
#
# Outputs (per dataset-countfile, under 05_score/transfer/):
#   - scores/<ID>__imrs_scores.tsv
#   - qc/<ID>__qc_summary.tsv
#   - qc/<ID>__top_contributors.tsv   (top genes by |w*z| per sample)
#
# Notes:
#   - Weights are loaded ONCE from anchors and reused for all scoring.
#   - Calibration datasets are scored identically; NO weight updates.
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
# CONFIG (PDF defaults)
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

weights_path <- file.path(project_root, "05_score", "anchors", "gene_weights.tsv")

# UPDATED: scoring design root (dataset-level)
design_root <- file.path(project_root, "00_metadata", "verified_metadata", "scoring")

min_ctrl_n       <- 3
gene_sd_floor    <- 0.10
score_sd_floor   <- 0.10
min_coverage     <- 0.80
strip_ens_ver    <- TRUE

# How many contributors to report PER sample
top_k_contrib <- 20

# Output base
out_base <- file.path(project_root, "05_score", "transfer")
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
  df <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(read_table(path, show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  stop("Failed to read design (delimiter/header issue): ", path)
}

# Find all gene_counts_clean.tsv under 03_counts/<DATASET>/
find_counts_files <- function(dataset_id) {
  root <- file.path(project_root, "03_counts", dataset_id)
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = "^gene_counts_clean\\.tsv$", recursive = TRUE, full.names = TRUE)
}

# Build a stable ID for outputs from dataset + relative path
make_counts_id <- function(dataset_id, counts_path) {
  root <- file.path(project_root, "03_counts", dataset_id)
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", gsub("\\\\", "/", counts_path))
  rel <- sub("/gene_counts_clean\\.tsv$", "", rel)
  rel <- ifelse(nchar(rel) == 0, "counts", rel)
  rel <- gsub("[^A-Za-z0-9_\\-]+", "_", rel)
  paste0(dataset_id, "__", rel)
}

# Read counts: first col gene_id, rest samples
read_counts_matrix <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (ncol(df) < 2) stop("Counts file has <2 columns: ", path)

  names(df)[1] <- "gene_id"
  gene_id <- as.character(df$gene_id)
  mat <- as.matrix(df[, -1, drop = FALSE])

  # Preserve numeric/integer values; validate integers later
  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- gene_id
  colnames(mat) <- trimws(colnames(mat))
  mat
}

# Strict integer check (PDF QC rule #1)
is_integerish_matrix <- function(mat) {
  all(is.finite(mat)) && all(abs(mat - round(mat)) < 1e-6)
}

# Load frozen weights (anchors only)
load_frozen_anchor_weights <- function(path) {
  if (!file.exists(path)) stop("Missing frozen anchor weights: ", path)
  w <- read_tsv(path, show_col_types = FALSE)

  # Accept either "gene" (spec) or "gene_id" (internal)
  gene_col <- NULL
  if ("gene_id" %in% names(w)) gene_col <- "gene_id"
  if (is.null(gene_col) && ("gene" %in% names(w))) gene_col <- "gene"
  if (is.null(gene_col)) {
    stop("gene_weights.tsv must contain a gene column named 'gene' (spec) or 'gene_id' (internal)")
  }

  # Accept either beta_meta (preferred) or weight (spec)
  if (!("beta_meta" %in% names(w)) && !("weight" %in% names(w))) {
    stop("gene_weights.tsv must contain 'beta_meta' (preferred) or 'weight'")
  }
  if (!("beta_meta" %in% names(w))) w$beta_meta <- w$weight

  # Optional flags (if present) can exclude genes
  if (!("heterogeneity_flag" %in% names(w))) w$heterogeneity_flag <- FALSE
  if (!("low_power_flag" %in% names(w)))      w$low_power_flag <- FALSE

  w <- w %>%
    transmute(
      gene_id = as.character(.data[[gene_col]]),
      beta_meta = suppressWarnings(as.numeric(beta_meta)),
      heterogeneity_flag = as.logical(heterogeneity_flag),
      low_power_flag = as.logical(low_power_flag)
    ) %>%
    filter(!is.na(beta_meta)) %>%
    filter(!(isTRUE(heterogeneity_flag) | isTRUE(low_power_flag)))

  if (strip_ens_ver) {
    w <- w %>% mutate(gene_id = strip_ensembl_version(gene_id))
  }

  w <- w %>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(beta_meta), n = 1, with_ties = FALSE) %>%
    ungroup()

  if (nrow(w) == 0) stop("No usable weights after filtering/cleaning: ", path)
  w
}

# -------------------------
# Core scoring function
# -------------------------
score_one_counts_file <- function(dataset_id, counts_path, design_df, weights_df) {
  counts_id <- make_counts_id(dataset_id, counts_path)

  counts <- read_counts_matrix(counts_path)

  if (strip_ens_ver) {
    rownames(counts) <- strip_ensembl_version(rownames(counts))
  }

  if (!("sample_id" %in% names(design_df))) stop("Design missing sample_id for dataset: ", dataset_id)
  if (!("condition_simple" %in% names(design_df))) {
    stop("Design missing condition_simple (CONTROL/DELIVERY) for dataset: ", dataset_id)
  }

  design <- design_df %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    ) %>%
    filter(sample_id %in% colnames(counts))

  if (nrow(design) == 0) {
    message("Skip (no overlapping samples): ", counts_id)
    return(invisible(NULL))
  }

  counts <- counts[, design$sample_id, drop = FALSE]

  if (!is_integerish_matrix(counts)) {
    message("FAIL (non-integer counts): ", counts_id, " | ", counts_path)
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

  sf <- estimateSizeFactorsForMatrix(counts_int)
  norm <- sweep(counts_int, 2, sf, "/")
  x <- log2(norm + 1)

  w_genes <- weights_df$gene_id
  overlap <- intersect(rownames(x), w_genes)
  coverage <- length(overlap) / length(w_genes)

  if (!is.finite(coverage) || coverage < min_coverage) {
    message("FAIL (coverage < ", min_coverage, "): ", counts_id,
            " | coverage=", sprintf("%.3f", coverage),
            " | overlap=", length(overlap), " / ", length(w_genes))
    qc <- tibble(
      id = counts_id,
      dataset_id = dataset_id,
      counts_path = counts_path,
      pass = FALSE,
      fail_reason = paste0("Coverage < ", min_coverage),
      coverage = coverage,
      n_weights = length(w_genes),
      n_overlap = length(overlap),
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
    imrs_z   = as.numeric(imrs_z[design$sample_id]),
    coverage = coverage
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
    n_weights = length(w_genes),
    n_overlap = length(overlap),
    coverage = coverage,
    n_genes_sd_floor = n_hit_floor,
    gene_sd_floor = gene_sd_floor,
    score_sd_floor = score_sd_floor,
    controls_imrs_z_mean = mean(imrs_z[ctrl_ids], na.rm = TRUE),
    controls_imrs_z_sd   = sd(imrs_z[ctrl_ids], na.rm = TRUE)
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

  message("OK: ", counts_id,
          " | n_samples=", nrow(design),
          " | n_ctrl=", length(ctrl_ids),
          " | coverage=", sprintf("%.3f", coverage),
          " | sd_floor_genes=", n_hit_floor)

  invisible(NULL)
}

# -------------------------
# MAIN
# -------------------------
if (!file.exists(weights_path)) stop("Frozen weights missing: ", weights_path)

weights_df <- load_frozen_anchor_weights(weights_path)
message("Loaded frozen anchor weights: n_genes=", nrow(weights_df),
        " | path=", weights_path)

counts_root <- file.path(project_root, "03_counts")
if (!dir.exists(counts_root)) stop("Missing 03_counts root: ", counts_root)

dataset_ids <- list.dirs(counts_root, recursive = FALSE, full.names = FALSE)
dataset_ids <- dataset_ids[grepl("^GSE", dataset_ids)]
if (length(dataset_ids) == 0) stop("No datasets found under: ", counts_root)

if (!dir.exists(design_root)) {
  stop("Missing design_root (run scoring-design builder first): ", design_root)
}

for (ds in dataset_ids) {
  # UPDATED: scoring design path
  design_path <- file.path(design_root, ds, paste0(ds, "_design.tsv"))

  if (!file.exists(design_path)) {
    message("Skip dataset (missing scoring design): ", ds, " | ", design_path)
    next
  }

  design_df <- read_table_robust(design_path)

  if (!("condition_simple" %in% names(design_df))) {
    message("Skip dataset (design missing condition_simple): ", ds)
    next
  }

  cfiles <- find_counts_files(ds)
  if (length(cfiles) == 0) {
    message("Skip dataset (no gene_counts_clean.tsv found): ", ds)
    next
  }

  message("\n==============================")
  message("DATASET: ", ds, " | count_files=", length(cfiles))
  message("==============================")

  for (cf in cfiles) {
    tryCatch(
      score_one_counts_file(ds, cf, design_df, weights_df),
      error = function(e) {
        message("ERROR scoring: ", ds, " | ", cf, " | ", conditionMessage(e))
      }
    )
  }
}

message("\nDONE: Step 08 scoring finished (STRICT PDF workflow; frozen anchor weights; no retraining).")