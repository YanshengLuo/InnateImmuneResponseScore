#!/usr/bin/env Rscript

# ============================================================
# IMRS Step 11 — Noise sensitivity / outlier removal rescore
#
# Purpose:
#   Re-score flagged datasets after removing extreme samples,
#   while preserving the exact Step 08 scoring workflow:
#     - same frozen anchor weights
#     - same DESeq2 size-factor normalization on all included samples
#     - same control-only gene z-scoring
#     - same gene SD floor
#     - same score SD floor
#     - no retraining
#
# Inputs:
#   1) <project_root>/05_score/failure_diagnosis/diagnosis/diagnosis_table.tsv
#      or
#      <project_root>/05_score/failure_diagnosis/low_imrs_datasets.tsv
#
#   2) Counts:
#      <project_root>/03_counts/<DATASET>/.../gene_counts_clean.tsv
#
#   3) Design:
#      <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#
#   4) Frozen weights:
#      <project_root>/05_score/anchors/gene_weights.tsv
#
# Outputs:
#   <project_root>/05_score/noise_rescore/
#       rescore_summary.tsv
#       rescore_sample_level.tsv
#       outlier_flags.tsv
#       rescored_scores/<ID>__rescored_scores.tsv
#       rescored_qc/<ID>__rescored_qc.tsv
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

diagnosis_file <- file.path(project_root, "05_score", "failure_diagnosis", "diagnosis", "diagnosis_table.tsv")
low_file       <- file.path(project_root, "05_score", "failure_diagnosis", "low_imrs_datasets.tsv")

weights_path <- file.path(project_root, "05_score", "anchors", "gene_weights.tsv")
design_root  <- file.path(project_root, "00_metadata", "verified_metadata", "scoring")

out_base       <- file.path(project_root, "05_score", "noise_rescore")
out_scores_dir <- file.path(out_base, "rescored_scores")
out_qc_dir     <- file.path(out_base, "rescored_qc")

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
dir.create(out_scores_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_qc_dir, showWarnings = FALSE, recursive = TRUE)

# exact Step 08 defaults
min_ctrl_n     <- 3
gene_sd_floor  <- 0.10
score_sd_floor <- 0.10
min_coverage   <- 0.80
strip_ens_ver  <- TRUE

# outlier settings
outlier_method <- "mad_by_group"   # "mad_by_group" or "mad_global"
mad_thresh <- 3.5
max_remove_per_group <- 1
max_remove_total <- 2

# selection of targets
selection_mode <- "noisy_only"   # "all_flagged", "noisy_only", "weak_and_noisy"

# -------------------------
# HELPERS
# -------------------------
strip_ensembl_version <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+$", "", x)
}

read_table_robust <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)

  df <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  df <- tryCatch(read_table(path, show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 2) return(df)

  stop("Failed to read table: ", path)
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

load_frozen_anchor_weights <- function(path) {
  if (!file.exists(path)) stop("Missing frozen anchor weights: ", path)
  w <- read_tsv(path, show_col_types = FALSE)

  gene_col <- NULL
  if ("gene_id" %in% names(w)) gene_col <- "gene_id"
  if (is.null(gene_col) && "gene" %in% names(w)) gene_col <- "gene"
  if (is.null(gene_col)) stop("gene_weights.tsv must contain gene_id or gene")

  if (!("beta_meta" %in% names(w)) && !("weight" %in% names(w))) {
    stop("gene_weights.tsv must contain beta_meta or weight")
  }
  if (!("beta_meta" %in% names(w))) w$beta_meta <- w$weight

  if (!("heterogeneity_flag" %in% names(w))) w$heterogeneity_flag <- FALSE
  if (!("low_power_flag" %in% names(w))) w$low_power_flag <- FALSE

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

  w %>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(beta_meta), n = 1, with_ties = FALSE) %>%
    ungroup()
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

mad_score <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1, na.rm = TRUE)
  if (!is.finite(madv) || madv <= 0) return(rep(0, length(x)))
  abs(x - med) / madv
}

find_counts_path_from_score_file <- function(score_file, project_root) {
  b <- basename(score_file)
  id <- sub("__imrs_scores\\.tsv$", "", b)
  dataset_id <- sub("__.*$", "", id)

  root <- file.path(project_root, "03_counts", dataset_id)
  if (!dir.exists(root)) return(NA_character_)

  candidates <- list.files(root, pattern = "^gene_counts_clean\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) == 0) return(NA_character_)

  make_counts_id <- function(dataset_id, counts_path) {
    root2 <- file.path(project_root, "03_counts", dataset_id)
    rel <- sub(paste0("^", gsub("\\\\", "/", root2), "/?"), "", gsub("\\\\", "/", counts_path))
    rel <- sub("/gene_counts_clean\\.tsv$", "", rel)
    rel <- ifelse(nchar(rel) == 0, "counts", rel)
    rel <- gsub("[^A-Za-z0-9_\\-]+", "_", rel)
    paste0(dataset_id, "__", rel)
  }

  ids <- vapply(candidates, function(x) make_counts_id(dataset_id, x), character(1))
  hit <- candidates[ids == id]
  if (length(hit) >= 1) return(hit[1])

  NA_character_
}

score_with_exact_step08_math <- function(dataset_id, counts_path, design_df, weights_df, keep_samples) {
  counts <- read_counts_matrix(counts_path)

  if (strip_ens_ver) {
    rownames(counts) <- strip_ensembl_version(rownames(counts))
  }

  design <- design_df %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    ) %>%
    filter(sample_id %in% colnames(counts)) %>%
    filter(sample_id %in% keep_samples)

  if (nrow(design) == 0) stop("No overlapping kept samples for: ", dataset_id)

  counts <- counts[, design$sample_id, drop = FALSE]

  if (!is_integerish_matrix(counts)) stop("Counts not integer for: ", dataset_id)

  counts_int <- round(counts)
  storage.mode(counts_int) <- "integer"

  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  if (length(ctrl_ids) < min_ctrl_n) stop("Controls < min_ctrl_n after filtering for: ", dataset_id)

  sf <- estimateSizeFactorsForMatrix(counts_int)
  norm <- sweep(counts_int, 2, sf, "/")
  x <- log2(norm + 1)

  w_genes <- weights_df$gene_id
  overlap <- intersect(rownames(x), w_genes)
  coverage <- length(overlap) / length(w_genes)

  if (!is.finite(coverage) || coverage < min_coverage) {
    stop("Coverage < min_coverage after filtering for: ", dataset_id)
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

  qc <- tibble(
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

  list(scores = out_scores, qc = qc)
}

flag_outliers <- function(score_df, method = "mad_by_group") {
  df <- score_df %>%
    mutate(
      sample_id = as.character(sample_id),
      condition_simple = toupper(as.character(condition_simple)),
      imrs_z = as.numeric(imrs_z)
    )

  if (method == "mad_global") {
    df <- df %>%
      mutate(
        robust_score = mad_score(imrs_z),
        outlier_flag = robust_score > mad_thresh
      )
  } else {
    df <- df %>%
      group_by(condition_simple) %>%
      mutate(
        robust_score = mad_score(imrs_z),
        outlier_flag = robust_score > mad_thresh
      ) %>%
      ungroup()
  }

  df <- df %>%
    group_by(condition_simple) %>%
    arrange(desc(outlier_flag), desc(robust_score), .by_group = TRUE) %>%
    mutate(remove_rank_group = ifelse(outlier_flag, row_number(), NA_integer_)) %>%
    ungroup() %>%
    arrange(desc(outlier_flag), desc(robust_score)) %>%
    mutate(remove_rank_total = ifelse(outlier_flag, row_number(), NA_integer_)) %>%
    mutate(
      remove_final = ifelse(
        isTRUE(outlier_flag) &
          !is.na(remove_rank_group) & remove_rank_group <= max_remove_per_group &
          !is.na(remove_rank_total) & remove_rank_total <= max_remove_total,
        TRUE, FALSE
      )
    )

  df
}

summarize_score_object <- function(score_df) {
  ctrl <- score_df %>% filter(condition_simple == "CONTROL")
  del  <- score_df %>% filter(condition_simple == "DELIVERY")

  tibble(
    n_control = nrow(ctrl),
    n_delivery = nrow(del),
    control_mean_imrs_z = safe_mean(ctrl$imrs_z),
    control_sd_imrs_z = safe_sd(ctrl$imrs_z),
    delivery_mean_imrs_z = safe_mean(del$imrs_z),
    delivery_sd_imrs_z = safe_sd(del$imrs_z),
    delta_imrs = safe_mean(del$imrs_z) - safe_mean(ctrl$imrs_z)
  )
}

# -------------------------
# LOAD TARGET TABLE
# -------------------------
if (file.exists(diagnosis_file)) {
  target_tbl <- read_tsv(diagnosis_file, show_col_types = FALSE, progress = FALSE)
  source_type <- "diagnosis"
} else if (file.exists(low_file)) {
  target_tbl <- read_tsv(low_file, show_col_types = FALSE, progress = FALSE)
  source_type <- "low"
} else {
  stop("Missing both diagnosis_table.tsv and low_imrs_datasets.tsv")
}

if (!all(c("dataset_id", "id", "score_file") %in% names(target_tbl))) {
  stop("Target table missing required columns: dataset_id, id, score_file")
}

if (source_type == "diagnosis" && "diagnosis" %in% names(target_tbl)) {
  if (selection_mode == "noisy_only") {
    target_tbl <- target_tbl %>% filter(diagnosis == "noisy_dataset")
  } else if (selection_mode == "weak_and_noisy") {
    target_tbl <- target_tbl %>% filter(diagnosis %in% c("noisy_dataset", "true_weak_biology"))
  }
}

if (nrow(target_tbl) == 0) {
  stop("No target rows selected for noise rescore.")
}

weights_df <- load_frozen_anchor_weights(weights_path)

summary_rows <- list()
sample_rows  <- list()
outlier_rows <- list()

# -------------------------
# MAIN LOOP
# -------------------------
for (i in seq_len(nrow(target_tbl))) {
  row <- target_tbl[i, , drop = FALSE]

  dataset_id <- as.character(row$dataset_id[[1]])
  id         <- as.character(row$id[[1]])
  score_file <- as.character(row$score_file[[1]])

  design_path <- file.path(design_root, dataset_id, paste0(dataset_id, "_design.tsv"))
  if (!file.exists(design_path)) {
    message("Skip (missing design): ", id)
    next
  }

  counts_path <- find_counts_path_from_score_file(score_file, project_root)
  if (is.na(counts_path) || !file.exists(counts_path)) {
    message("Skip (could not resolve counts path): ", id)
    next
  }

  design_df <- read_table_robust(design_path)
  score_df  <- read_tsv(score_file, show_col_types = FALSE, progress = FALSE)

  if (!all(c("sample_id", "condition_simple", "imrs_z") %in% names(score_df))) {
    message("Skip (score file missing columns): ", id)
    next
  }

  score_df <- score_df %>%
    mutate(
      sample_id = as.character(sample_id),
      condition_simple = toupper(as.character(condition_simple)),
      imrs_z = as.numeric(imrs_z)
    )

  before_sum <- summarize_score_object(score_df)

  flagged <- flag_outliers(score_df, method = outlier_method)
  keep_samples <- flagged %>%
    filter(!remove_final) %>%
    pull(sample_id)

  removed_samples <- flagged %>%
    filter(remove_final)

  outlier_rows[[length(outlier_rows) + 1]] <- flagged %>%
    mutate(
      dataset_id = dataset_id,
      id = id,
      counts_path = counts_path,
      score_file = score_file
    )

  if (nrow(removed_samples) == 0) {
    after_scores <- score_df
    after_qc <- tibble(
      dataset_id = dataset_id,
      counts_path = counts_path,
      pass = TRUE,
      note = "No outliers removed"
    )
    rescore_status <- "no_removal"
    remove_n <- 0
  } else {
    rescored <- tryCatch(
      score_with_exact_step08_math(
        dataset_id = dataset_id,
        counts_path = counts_path,
        design_df = design_df,
        weights_df = weights_df,
        keep_samples = keep_samples
      ),
      error = function(e) e
    )

    if (inherits(rescored, "error")) {
      message("Rescore failed: ", id, " | ", conditionMessage(rescored))

      summary_row <- bind_cols(
        tibble(
          dataset_id = dataset_id,
          id = id,
          counts_path = counts_path,
          score_file = score_file,
          design_file = design_path,
          outlier_method = outlier_method,
          mad_thresh = mad_thresh,
          remove_n = nrow(removed_samples),
          removed_samples = paste(removed_samples$sample_id, collapse = ";"),
          rescore_status = "rescore_failed",
          fail_reason = conditionMessage(rescored)
        ),
        before_sum %>%
          rename_with(~ paste0("before_", .x)),
        tibble(
          after_n_control = NA_integer_,
          after_n_delivery = NA_integer_,
          after_control_mean_imrs_z = NA_real_,
          after_control_sd_imrs_z = NA_real_,
          after_delivery_mean_imrs_z = NA_real_,
          after_delivery_sd_imrs_z = NA_real_,
          after_delta_imrs = NA_real_,
          delta_imrs_change = NA_real_,
          delivery_sd_change = NA_real_
        )
      )

      summary_rows[[length(summary_rows) + 1]] <- summary_row
      next
    }

    after_scores <- rescored$scores
    after_qc     <- rescored$qc
    rescore_status <- "rescored"
    remove_n <- nrow(removed_samples)

    write_tsv(
      after_scores,
      file.path(out_scores_dir, paste0(id, "__rescored_scores.tsv"))
    )
    write_tsv(
      after_qc,
      file.path(out_qc_dir, paste0(id, "__rescored_qc.tsv"))
    )
  }

  after_sum <- summarize_score_object(after_scores)

  summary_row <- bind_cols(
    tibble(
      dataset_id = dataset_id,
      id = id,
      counts_path = counts_path,
      score_file = score_file,
      design_file = design_path,
      outlier_method = outlier_method,
      mad_thresh = mad_thresh,
      remove_n = remove_n,
      removed_samples = ifelse(remove_n > 0, paste(removed_samples$sample_id, collapse = ";"), NA_character_),
      rescore_status = rescore_status,
      fail_reason = NA_character_
    ),
    before_sum %>%
      rename_with(~ paste0("before_", .x)),
    after_sum %>%
      rename_with(~ paste0("after_", .x))
  ) %>%
    mutate(
      delta_imrs_change = after_delta_imrs - before_delta_imrs,
      delivery_sd_change = after_delivery_sd_imrs_z - before_delivery_sd_imrs_z
    )

  summary_rows[[length(summary_rows) + 1]] <- summary_row

  sample_rows[[length(sample_rows) + 1]] <- bind_rows(
    score_df %>%
      mutate(
        dataset_id = dataset_id,
        id = id,
        version = "before"
      ),
    after_scores %>%
      mutate(
        dataset_id = dataset_id,
        id = id,
        version = "after"
      )
  )
}

summary_tbl <- bind_rows(summary_rows)
sample_tbl  <- bind_rows(sample_rows)
outlier_tbl <- bind_rows(outlier_rows)

write_tsv(summary_tbl, file.path(out_base, "rescore_summary.tsv"))
write_tsv(sample_tbl,  file.path(out_base, "rescore_sample_level.tsv"))
write_tsv(outlier_tbl, file.path(out_base, "outlier_flags.tsv"))

message("Saved:")
message("  ", file.path(out_base, "rescore_summary.tsv"))
message("  ", file.path(out_base, "rescore_sample_level.tsv"))
message("  ", file.path(out_base, "outlier_flags.tsv"))
message("  ", out_scores_dir)
message("  ", out_qc_dir)
message("Done.")