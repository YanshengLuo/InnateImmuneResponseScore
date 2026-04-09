#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ============================================================
# 04_compute_drs_components.R
#
# Two-pass prototype:
#   Pass 1: compute raw per-dataset metrics
#   Pass 2: calibrate component scores from empirical distributions
#
# Updated:
#   - replaces min-component failure mode assignment
#   - uses rule-based decision tree for prototype diagnosis
#   - adds prototype_bio_label
#   - optionally merges a manual truth table for comparison
#
# Outputs:
#   intermediate/04_drs_raw_metrics.tsv
#   intermediate/04_drs_component_metrics.tsv
#   intermediate/04_drs_per_sample_diagnostics.tsv
#   intermediate/04_drs_calibration_bounds.tsv
# ============================================================

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
framework_root <- if (length(args) >= 1) args[1] else
  "D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/Dataset_Level_Confidence_score_prototype"

inputs_dir <- file.path(framework_root, "inputs")
intermediate_dir <- file.path(framework_root, "intermediate")
logs_dir <- file.path(framework_root, "logs")

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(inputs_dir, "01_dataset_manifest.tsv")
sample_file   <- file.path(inputs_dir, "02_sample_metadata.tsv")
index_file    <- file.path(inputs_dir, "03_expression_matrix_index.tsv")

raw_metrics_out <- file.path(intermediate_dir, "04_drs_raw_metrics.tsv")
component_out   <- file.path(intermediate_dir, "04_drs_component_metrics.tsv")
sample_out      <- file.path(intermediate_dir, "04_drs_per_sample_diagnostics.tsv")
bounds_out      <- file.path(intermediate_dir, "04_drs_calibration_bounds.tsv")
log_file        <- file.path(logs_dir, "04_compute_drs_components.log")

# OPTIONAL: manually set this to your truth table path
truth_table_file <- file.path(framework_root, "inputs", "imrs_prototype_truth_table.tsv")

if (file.exists(log_file)) file.remove(log_file)

# prototype defaults
prior_count <- 0.5
min_detected_samples <- 2
top_k_genes <- 2000
corr_method <- "spearman"
pca_var_cutoff <- 0.80
pca_max_k <- 10
outlier_z_thresh <- 3
min_group_n_for_geometry <- 2

# context defaults
ctx_n0 <- 2
ctx_n1 <- 6
ctx_v0 <- 0.20
ctx_v1 <- 0.80

# softened outlier penalty
outlier_alpha <- 2
outlier_beta <- 0.5

# empirical calibration quantiles
q_low <- 0.10
q_high <- 0.90

# diagnosis thresholds
diag_outlier_fraction <- 0.20
diag_outlier_severity <- 1.0
diag_low_context_n <- 3
diag_low_coh <- 0.35
diag_low_sep <- 0.35
diag_mid_sep <- 0.45
diag_mid_coh <- 0.45

# -------------------------
# LOGGING
# -------------------------
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

# -------------------------
# HELPERS
# -------------------------
clip01 <- function(x) pmax(0, pmin(1, x))

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
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

safe_mad <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mad(x, constant = 1)
}

robust_z <- function(x) {
  med <- safe_median(x)
  madv <- safe_mad(x)
  if (!is.finite(madv) || madv <= 0) return(rep(0, length(x)))
  (x - med) / (1.4826 * madv)
}

norm_entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) <= 1) return(0)
  p <- x / sum(x)
  h <- -sum(p * log(p))
  hmax <- log(length(p))
  if (!is.finite(hmax) || hmax <= 0) return(0)
  h / hmax
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

make_bounds <- function(x, low = 0.10, high = 0.90, fallback_low = NA_real_, fallback_high = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x) < 3) {
    return(list(low = fallback_low, high = fallback_high))
  }
  lo <- safe_quantile(x, low)
  hi <- safe_quantile(x, high)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    lo <- min(x, na.rm = TRUE)
    hi <- max(x, na.rm = TRUE)
  }
  list(low = lo, high = hi)
}

scale_positive_metric <- function(x, lo, hi) {
  if (!is.finite(x) || !is.finite(lo) || !is.finite(hi) || hi <= lo) return(NA_real_)
  clip01((x - lo) / (hi - lo))
}

scale_negative_metric <- function(x, lo, hi) {
  if (!is.finite(x) || !is.finite(lo) || !is.finite(hi) || hi <= lo) return(NA_real_)
  clip01((hi - x) / (hi - lo))
}

read_expr_matrix <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  stopifnot(ncol(df) >= 2)
  names(df)[1] <- "gene_id"

  gene_id <- as.character(df$gene_id)
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_id
  colnames(mat) <- colnames(df)[-1]
  mat
}

logcpm_transform <- function(counts, prior_count = 0.5) {
  lib_size <- colSums(counts, na.rm = TRUE)
  lib_size[lib_size <= 0 | !is.finite(lib_size)] <- NA_real_
  cpm <- sweep(counts, 2, lib_size / 1e6, "/")
  log2(cpm + prior_count)
}

filter_and_select_genes <- function(expr_log, k = 2000, min_samples = 2) {
  keep_detected <- rowSums(is.finite(expr_log) & expr_log > 0, na.rm = TRUE) >= min_samples
  x <- expr_log[keep_detected, , drop = FALSE]

  if (nrow(x) == 0) return(x)

  madv <- apply(x, 1, safe_mad)
  madv[!is.finite(madv)] <- -Inf
  ord <- order(madv, decreasing = TRUE)

  k_use <- min(k, nrow(x))
  x[ord[seq_len(k_use)], , drop = FALSE]
}

standardize_genewise <- function(expr_log) {
  mu <- apply(expr_log, 1, safe_median)
  madv <- apply(expr_log, 1, safe_mad)
  s <- 1.4826 * madv
  s[!is.finite(s) | s <= 0] <- 1

  z <- sweep(expr_log, 1, mu, "-")
  z <- sweep(z, 1, s, "/")
  z
}

choose_k_pcs <- function(var_explained, cutoff = 0.80, max_k = 10) {
  if (length(var_explained) == 0) return(1L)
  cumv <- cumsum(var_explained)
  k80 <- which(cumv >= cutoff)[1]
  if (is.na(k80)) k80 <- length(var_explained)
  as.integer(max(1, min(k80, max_k, length(var_explained))))
}

compute_pca_scores <- function(z_gene_sample, cutoff = 0.80, max_k = 10) {
  x <- t(z_gene_sample)
  x <- x[, apply(x, 2, function(v) any(is.finite(v))), drop = FALSE]

  if (nrow(x) < 2 || ncol(x) < 2) {
    return(list(
      scores = matrix(0, nrow = nrow(x), ncol = 1, dimnames = list(rownames(x), "PC1")),
      var_explained = 1
    ))
  }

  pca <- prcomp(x, center = TRUE, scale. = FALSE)
  ve <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  k <- choose_k_pcs(ve, cutoff = cutoff, max_k = max_k)
  scores <- pca$x[, seq_len(k), drop = FALSE]
  list(scores = scores, var_explained = ve)
}

compute_sample_correlation <- function(z_gene_sample, method = "spearman") {
  suppressWarnings(cor(z_gene_sample, method = method, use = "pairwise.complete.obs"))
}

pairwise_group_median_corr <- function(cor_mat, groups) {
  ug <- unique(groups)
  vals <- map_dbl(ug, function(g) {
    idx <- which(groups == g)
    if (length(idx) < 2) return(NA_real_)
    cm <- cor_mat[idx, idx, drop = FALSE]
    safe_median(cm[upper.tri(cm, diag = FALSE)])
  })
  tibble(group = ug, median_corr = vals)
}

within_group_spread <- function(pca_scores, groups) {
  ug <- unique(groups)
  vals <- map_dbl(ug, function(g) {
    idx <- which(groups == g)
    if (length(idx) < 2) return(NA_real_)
    pts <- pca_scores[idx, , drop = FALSE]
    ctr <- colMeans(pts)
    d <- sqrt(rowSums((sweep(pts, 2, ctr, "-")) ^ 2))
    safe_median(d)
  })
  tibble(group = ug, median_spread = vals)
}

compute_sample_coherence_raw <- function(z_gene_sample, sample_meta) {
  cor_mat <- compute_sample_correlation(z_gene_sample, method = corr_method)
  pca_obj <- compute_pca_scores(z_gene_sample, cutoff = pca_var_cutoff, max_k = pca_max_k)

  group_corr <- pairwise_group_median_corr(cor_mat, sample_meta$group)
  group_spread <- within_group_spread(pca_obj$scores, sample_meta$group)

  rwithin <- safe_median(group_corr$median_corr)
  spreadwithin <- safe_median(group_spread$median_spread)

  list(
    raw = tibble(
      coherence_rwithin = rwithin,
      coherence_spreadwithin = spreadwithin
    ),
    per_sample = tibble(
      sample_id = rownames(pca_obj$scores),
      pca_pc1 = pca_obj$scores[, 1],
      pca_pc2 = if (ncol(pca_obj$scores) >= 2) pca_obj$scores[, 2] else 0
    ),
    objects = list(cor_mat = cor_mat, pca_scores = pca_obj$scores)
  )
}

compute_outlier_burden <- function(pca_scores, cor_mat, sample_meta) {
  ctr <- colMeans(pca_scores)
  dist_ctr <- sqrt(rowSums((sweep(pca_scores, 2, ctr, "-")) ^ 2))
  z_dist <- robust_z(dist_ctr)

  med_corr <- map_dbl(seq_len(nrow(cor_mat)), function(i) {
    others <- setdiff(seq_len(nrow(cor_mat)), i)
    safe_median(cor_mat[i, others])
  })
  z_conn <- robust_z(med_corr)

  outlier_flag <- (z_dist > outlier_z_thresh) | (z_conn < -outlier_z_thresh)

  outlier_severity <- pmax(0, z_dist - outlier_z_thresh, -z_conn - outlier_z_thresh)
  f_out <- mean(outlier_flag, na.rm = TRUE)
  sev <- safe_median(outlier_severity[outlier_flag])
  if (!is.finite(sev)) sev <- 0

  c_out <- exp(-outlier_alpha * f_out) * exp(-outlier_beta * sev)

  per_sample <- tibble(
    sample_id = rownames(pca_scores),
    outlier_dist = dist_ctr,
    outlier_z_dist = z_dist,
    median_corr_to_others = med_corr,
    outlier_z_connectivity = z_conn,
    outlier_flag = outlier_flag
  ) %>%
    left_join(sample_meta %>% select(sample_id, group), by = "sample_id")

  list(
    component = tibble(
      Cout = c_out,
      outlier_fraction = f_out,
      outlier_severity = sev
    ),
    per_sample = per_sample
  )
}

compute_silhouette_simple <- function(dist_mat, groups) {
  n <- nrow(dist_mat)
  out <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    gi <- groups[i]
    same <- which(groups == gi & seq_len(n) != i)
    other_groups <- setdiff(unique(groups), gi)

    if (length(same) == 0 || length(other_groups) == 0) {
      out[i] <- NA_real_
      next
    }

    a_i <- mean(dist_mat[i, same], na.rm = TRUE)

    b_i <- min(map_dbl(other_groups, function(g) {
      idx <- which(groups == g)
      mean(dist_mat[i, idx], na.rm = TRUE)
    }), na.rm = TRUE)

    out[i] <- (b_i - a_i) / max(a_i, b_i)
  }

  out
}

compute_group_separability_raw <- function(pca_scores, cor_mat, sample_meta) {
  groups <- sample_meta$group
  ug <- unique(groups)

  centroids <- map(ug, function(g) {
    idx <- which(groups == g)
    colMeans(pca_scores[idx, , drop = FALSE])
  })
  names(centroids) <- ug

  within_spreads <- map_dbl(ug, function(g) {
    idx <- which(groups == g)
    if (length(idx) < 2) return(NA_real_)
    pts <- pca_scores[idx, , drop = FALSE]
    ctr <- colMeans(pts)
    d <- sqrt(rowSums((sweep(pts, 2, ctr, "-")) ^ 2))
    safe_median(d)
  })

  pair_dists <- c()
  if (length(ug) >= 2) {
    combs <- combn(ug, 2, simplify = FALSE)
    pair_dists <- map_dbl(combs, function(cc) {
      sqrt(sum((centroids[[cc[1]]] - centroids[[cc[2]]]) ^ 2))
    })
  }

  between_dist <- safe_mean(pair_dists)
  within_dist <- safe_median(within_spreads)
  sep_ratio <- between_dist / ifelse(!is.finite(within_dist) || within_dist <= 0, 1e-6, within_dist)

  dist_mat <- 1 - cor_mat
  sil <- compute_silhouette_simple(dist_mat, groups)
  sil_mean <- safe_mean(sil)

  list(
    raw = tibble(
      separation_ratio = sep_ratio,
      silhouette_mean = sil_mean,
      between_centroid_distance = between_dist,
      within_group_spread = within_dist
    ),
    per_sample = tibble(
      sample_id = rownames(cor_mat),
      silhouette = sil
    )
  )
}

compute_cramers_v <- function(group, batch) {
  tab <- table(group, batch)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)

  suppressWarnings({
    cs <- chisq.test(tab, correct = FALSE)
  })
  n <- sum(tab)
  if (!is.finite(n) || n <= 0) return(NA_real_)

  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  if (k <= 0) return(NA_real_)

  sqrt(as.numeric(cs$statistic) / (n * k))
}

compute_context_sufficiency <- function(sample_meta, manifest_row) {
  grp_tab <- table(sample_meta$group)
  n_min <- if (length(grp_tab) == 0) 0 else min(grp_tab)

  c_rep <- clip01((n_min - ctx_n0) / (ctx_n1 - ctx_n0))
  c_bal <- norm_entropy(as.numeric(grp_tab))

  batch_vec <- sample_meta$batch
  batch_vec[batch_vec == ""] <- NA
  if (all(is.na(batch_vec))) {
    c_batch <- 1
    batch_v <- NA_real_
  } else {
    batch_v <- compute_cramers_v(sample_meta$group, batch_vec)
    if (!is.finite(batch_v)) {
      c_batch <- 1
    } else {
      c_batch <- 1 - clip01((batch_v - ctx_v0) / (ctx_v1 - ctx_v0))
    }
  }

  required_fields <- c(
    !is.na(manifest_row$organism) && manifest_row$organism != "",
    !is.na(manifest_row$assay_type) && manifest_row$assay_type != "",
    nrow(sample_meta) > 0,
    all(!is.na(sample_meta$group) & sample_meta$group != "")
  )
  c_meta <- mean(required_fields)

  c_ctx <- 0.35 * c_rep + 0.20 * c_bal + 0.25 * c_batch + 0.20 * c_meta

  tibble(
    Cctx = c_ctx,
    context_replicate_score = c_rep,
    context_balance_score = c_bal,
    context_batch_score = c_batch,
    context_metadata_score = c_meta,
    min_group_n = n_min,
    n_groups = length(grp_tab),
    batch_cramers_v = batch_v
  )
}

# -------------------------
# INTERNAL COMPONENT DIAGNOSIS
# -------------------------
assign_primary_failure_mode <- function(Ccoh, Cout, Csep, Cctx,
                                        outlier_fraction, outlier_severity,
                                        coherence_rwithin, coherence_spreadwithin,
                                        separation_ratio, silhouette_mean,
                                        min_group_n) {

  if (is.finite(outlier_fraction) &&
      is.finite(outlier_severity) &&
      outlier_fraction >= diag_outlier_fraction &&
      outlier_severity >= diag_outlier_severity) {
    return("Outlier")
  }

  if (is.finite(min_group_n) && min_group_n < diag_low_context_n) {
    return("Context")
  }

  if (is.finite(Ccoh) && is.finite(Csep) &&
      Ccoh < diag_low_coh && Csep < diag_mid_sep) {
    return("Coherence")
  }

  if (is.finite(Csep) && Csep < diag_low_sep &&
      is.finite(Ccoh) && Ccoh >= diag_low_coh) {
    return("Separation")
  }

  if (is.finite(Ccoh) && is.finite(Csep) &&
      Ccoh < diag_mid_coh && Csep < diag_low_sep) {
    return("Coherence")
  }

  vals <- c(
    Context    = unname(as.numeric(Cctx)[1]),
    Coherence  = unname(as.numeric(Ccoh)[1]),
    Separation = unname(as.numeric(Csep)[1]),
    Outlier    = unname(as.numeric(Cout)[1])
  )
  names(vals)[which.min(vals)][1]
}

# -------------------------
# FINAL BIOLOGICAL LABEL
# -------------------------
assign_prototype_bio_label <- function(Ccoh, Cout, Csep, Cctx,
                                       outlier_fraction, outlier_severity,
                                       coherence_rwithin, coherence_spreadwithin,
                                       separation_ratio, silhouette_mean,
                                       min_group_n,
                                       delta_imrs = NA_real_,
                                       delivery_sd_imrs_z = NA_real_) {

  # 1. technical / outlier issue:
  # allow either explicit outlier burden OR globally unstable structure
  if ((is.finite(outlier_fraction) && outlier_fraction >= 0.20 &&
       is.finite(outlier_severity) && outlier_severity >= 0.8) ||
      ((is.finite(Ccoh) && Ccoh < 0.25) &&
       (is.finite(Csep) && Csep < 0.25))) {
    return("Technical / outlier issue")
  }

  # 2. strong biology:
  # primary strong rule made stricter to avoid overcalling
  if (is.finite(Csep) && Csep >= 0.75 &&
      is.finite(Ccoh) && Ccoh >= 0.45 &&
      is.finite(Cout) && Cout >= 0.70) {
    return("Strong biology")
  }

  # 3. alternate strong rule:
  # rescue genuinely strong datasets with excellent separation
  if (is.finite(Csep) && Csep >= 0.85 &&
      is.finite(Cout) && Cout >= 0.70) {
    return("Strong biology")
  }

  # 4. heterogeneous biology:
  # loosened to catch distributed variability instead of miscalling weak
  if (is.finite(Ccoh) && Ccoh < 0.55 &&
      is.finite(Csep) && Csep >= 0.20 &&
      is.finite(Cout) && Cout >= 0.55) {
    return("Heterogeneous biology")
  }

  # 5. weak biology:
  # weak separation without strong evidence of technical failure
  if (is.finite(Csep) && Csep < 0.45 &&
      is.finite(Cout) && Cout >= 0.50) {
    return("Weak biology")
  }

  # fallbacks
  if (is.finite(Csep) && Csep >= 0.70) return("Strong biology")
  if (is.finite(Ccoh) && Ccoh < 0.55) return("Heterogeneous biology")
  if (is.finite(Cout) && Cout < 0.50) return("Technical / outlier issue")

  "Weak biology"
}

compare_to_truth <- function(prototype_label, truth_label) {
  case_when(
    is.na(truth_label) | truth_label == "" ~ NA_character_,
    is.na(prototype_label) | prototype_label == "" ~ "missing_prototype_label",
    prototype_label == truth_label ~ "match",
    prototype_label == "Technical / outlier issue" & truth_label != "Technical / outlier issue" ~ "overcalled_outlier",
    prototype_label == "Heterogeneous biology" & truth_label == "Weak biology" ~ "heterogeneous_vs_weak",
    prototype_label == "Weak biology" & truth_label == "Heterogeneous biology" ~ "weak_vs_heterogeneous",
    prototype_label == "Strong biology" & truth_label != "Strong biology" ~ "overcalled_strong",
    TRUE ~ paste0("mismatch:", prototype_label, "_vs_", truth_label)
  )
}

# -------------------------
# INPUT CHECK
# -------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(sample_file))
stopifnot(file.exists(index_file))

manifest <- read_tsv(manifest_file, show_col_types = FALSE, progress = FALSE)
sample_meta_all <- read_tsv(sample_file, show_col_types = FALSE, progress = FALSE)
index_df <- read_tsv(index_file, show_col_types = FALSE, progress = FALSE)

required_manifest <- c("dataset_id", "dataset_class", "organism", "assay_type", "treatment_type", "include_in_analysis")
required_sample <- c("dataset_id", "sample_id", "group", "replicate")
required_index <- c("dataset_id", "expression_matrix_file", "value_type", "gene_id_type", "n_genes", "n_samples")

stopifnot(all(required_manifest %in% names(manifest)))
stopifnot(all(required_sample %in% names(sample_meta_all)))
stopifnot(all(required_index %in% names(index_df)))

# -------------------------
# PASS 1: RAW METRICS
# -------------------------
raw_rows <- list()
sample_rows <- list()

datasets <- manifest %>%
  filter(toupper(as.character(include_in_analysis)) == "TRUE") %>%
  pull(dataset_id) %>%
  unique()

log_message("Datasets to process: ", length(datasets))

for (dataset_id in datasets) {
  log_message("Processing dataset: ", dataset_id)

  man_row <- manifest %>% filter(dataset_id == !!dataset_id) %>% slice(1)
  idx_row <- index_df %>% filter(dataset_id == !!dataset_id) %>% slice(1)
  sm <- sample_meta_all %>% filter(dataset_id == !!dataset_id)

  if (nrow(idx_row) == 0 || nrow(sm) == 0) {
    log_message("Skipping ", dataset_id, " because index/sample metadata missing.")
    next
  }

  expr_path <- file.path(framework_root, idx_row$expression_matrix_file[[1]])
  if (!file.exists(expr_path)) {
    log_message("Missing expression matrix for ", dataset_id, ": ", expr_path)
    next
  }

  if (idx_row$value_type[[1]] != "raw_counts") {
    log_message("Skipping ", dataset_id, " because prototype expects raw_counts.")
    next
  }

  expr_counts <- tryCatch(read_expr_matrix(expr_path), error = function(e) e)
  if (inherits(expr_counts, "error")) {
    log_message("Read failed for ", dataset_id, ": ", conditionMessage(expr_counts))
    next
  }

  keep_samples <- intersect(sm$sample_id, colnames(expr_counts))
  if (length(keep_samples) < 4) {
    log_message("Skipping ", dataset_id, " because too few overlapping samples.")
    next
  }

  sm <- sm %>%
    filter(sample_id %in% keep_samples) %>%
    mutate(group = as.character(group)) %>%
    arrange(match(sample_id, keep_samples))

  expr_counts <- expr_counts[, sm$sample_id, drop = FALSE]

  grp_tab <- table(sm$group)
  if (length(grp_tab) < 2 || any(grp_tab < min_group_n_for_geometry)) {
    log_message("Skipping ", dataset_id, " because groups insufficient for geometry.")
    next
  }

  expr_log <- logcpm_transform(expr_counts, prior_count = prior_count)
  expr_sel <- filter_and_select_genes(expr_log, k = top_k_genes, min_samples = min_detected_samples)

  if (nrow(expr_sel) < 50) {
    log_message("Skipping ", dataset_id, " because too few informative genes after filtering.")
    next
  }

  z_gene_sample <- standardize_genewise(expr_sel)

  coh <- compute_sample_coherence_raw(z_gene_sample, sm)
  out <- compute_outlier_burden(coh$objects$pca_scores, coh$objects$cor_mat, sm)
  sep <- compute_group_separability_raw(coh$objects$pca_scores, coh$objects$cor_mat, sm)
  ctx <- compute_context_sufficiency(sm, man_row)

  raw_row <- bind_cols(
    tibble(
      dataset_id = dataset_id,
      dataset_class = man_row$dataset_class[[1]],
      organism = man_row$organism[[1]],
      assay_type = man_row$assay_type[[1]],
      tissue = if ("tissue" %in% names(man_row)) man_row$tissue[[1]] else NA_character_,
      treatment_type = man_row$treatment_type[[1]],
      n_samples_used = ncol(expr_counts),
      n_genes_selected = nrow(expr_sel)
    ),
    coh$raw,
    out$component,
    sep$raw,
    ctx
  )

  raw_rows[[length(raw_rows) + 1]] <- raw_row

  per_sample <- sm %>%
    select(dataset_id, sample_id, group, replicate, batch, time_hr = any_of("time_hr")) %>%
    left_join(coh$per_sample, by = "sample_id") %>%
    left_join(out$per_sample %>% select(-group), by = "sample_id") %>%
    left_join(sep$per_sample, by = "sample_id")

  sample_rows[[length(sample_rows) + 1]] <- per_sample

  log_message(
    "Finished raw metrics for ", dataset_id,
    " | rwithin=", round(raw_row$coherence_rwithin, 3),
    " spread=", round(raw_row$coherence_spreadwithin, 3),
    " out_frac=", round(raw_row$outlier_fraction, 3),
    " sep_ratio=", round(raw_row$separation_ratio, 3),
    " sil=", round(raw_row$silhouette_mean, 3)
  )
}

raw_tbl <- bind_rows(raw_rows) %>%
  arrange(dataset_class, dataset_id)

sample_tbl <- bind_rows(sample_rows) %>%
  arrange(dataset_id, group, replicate, sample_id)

if (nrow(raw_tbl) == 0) {
  stop("No datasets were successfully processed.")
}

write_tsv(raw_tbl, raw_metrics_out, na = "")
write_tsv(sample_tbl, sample_out, na = "")

# -------------------------
# PASS 2: EMPIRICAL CALIBRATION
# -------------------------
b_coh_cor    <- make_bounds(raw_tbl$coherence_rwithin, low = q_low, high = q_high)
b_coh_spread <- make_bounds(raw_tbl$coherence_spreadwithin, low = q_low, high = q_high)
b_sep_ratio  <- make_bounds(log1p(raw_tbl$separation_ratio), low = q_low, high = q_high)
b_sep_sil    <- make_bounds(raw_tbl$silhouette_mean, low = q_low, high = q_high)

bounds_tbl <- tibble(
  metric = c("coherence_rwithin", "coherence_spreadwithin", "log1p_separation_ratio", "silhouette_mean"),
  low = c(b_coh_cor$low, b_coh_spread$low, b_sep_ratio$low, b_sep_sil$low),
  high = c(b_coh_cor$high, b_coh_spread$high, b_sep_ratio$high, b_sep_sil$high),
  q_low = q_low,
  q_high = q_high
)

component_tbl <- raw_tbl %>%
  mutate(
    coherence_corr_score = map_dbl(coherence_rwithin, ~ scale_positive_metric(.x, b_coh_cor$low, b_coh_cor$high)),
    coherence_spread_score = map_dbl(coherence_spreadwithin, ~ scale_negative_metric(.x, b_coh_spread$low, b_coh_spread$high)),
    Ccoh = (coherence_corr_score + coherence_spread_score) / 2,

    separation_ratio_score = map_dbl(log1p(separation_ratio), ~ scale_positive_metric(.x, b_sep_ratio$low, b_sep_ratio$high)),
    separation_silhouette_score = map_dbl(silhouette_mean, ~ scale_positive_metric(.x, b_sep_sil$low, b_sep_sil$high)),
    Csep = (separation_ratio_score + separation_silhouette_score) / 2
  ) %>%
  mutate(
    primary_failure_mode = pmap_chr(
      list(Ccoh, Cout, Csep, Cctx,
           outlier_fraction, outlier_severity,
           coherence_rwithin, coherence_spreadwithin,
           separation_ratio, silhouette_mean,
           min_group_n),
      assign_primary_failure_mode
    ),
    prototype_bio_label = pmap_chr(
      list(Ccoh, Cout, Csep, Cctx,
           outlier_fraction, outlier_severity,
           coherence_rwithin, coherence_spreadwithin,
           separation_ratio, silhouette_mean,
           min_group_n),
      assign_prototype_bio_label
    )
  ) %>%
  arrange(dataset_class, dataset_id)

# -------------------------
# OPTIONAL TRUTH TABLE MERGE
# -------------------------
if (file.exists(truth_table_file)) {
  truth_tbl <- read_tsv(truth_table_file, show_col_types = FALSE, progress = FALSE)

  truth_label_col <- NULL
  candidates <- c("manual_label", "truth_label", "manual_class", "label")
  hit <- candidates[candidates %in% names(truth_tbl)]
  if (length(hit) > 0) truth_label_col <- hit[1]

  if (!is.null(truth_label_col) && !is.na(truth_label_col)) {
    truth_tbl2 <- truth_tbl %>%
      transmute(
        dataset_id = as.character(dataset_id),
        truth_table_label = as.character(.data[[truth_label_col]])
      )

    component_tbl <- component_tbl %>%
      left_join(truth_tbl2, by = "dataset_id") %>%
      mutate(
        prototype_label_match = compare_to_truth(prototype_bio_label, truth_table_label)
      )
  } else {
    log_message("Truth table found, but no recognized label column.")
  }
} else {
  log_message("No truth table file found. Proceeding without truth-table comparison.")
}

write_tsv(bounds_tbl, bounds_out, na = "")
write_tsv(component_tbl, component_out, na = "")

log_message("Saved raw metrics: ", raw_metrics_out)
log_message("Saved calibration bounds: ", bounds_out)
log_message("Saved component metrics: ", component_out)
log_message("Saved per-sample diagnostics: ", sample_out)
log_message("Datasets scored: ", nrow(component_tbl))
log_message("Done.")