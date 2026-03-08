#!/usr/bin/env Rscript
# ============================================================
# IMRS threshold sensitivity summary
#
# Reads threshold snapshots from:
#   <project_root>/05_score/sensitivity_runs/thr_*/...
#
# Writes:
#   <project_root>/05_score/sensitivity_runs/threshold_sensitivity_summary.tsv
#   <project_root>/05_score/sensitivity_runs/threshold_sensitivity_by_dataset.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

sens_root <- file.path(project_root, "05_score", "sensitivity_runs")
if (!dir.exists(sens_root)) {
  stop("Missing sensitivity_runs folder: ", sens_root)
}

threshold_dirs <- list.dirs(sens_root, recursive = FALSE, full.names = TRUE)
threshold_dirs <- threshold_dirs[grepl("thr_", basename(threshold_dirs))]

if (length(threshold_dirs) == 0) {
  stop("No threshold snapshot folders found under: ", sens_root)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

read_one_run <- function(run_dir) {
  thr_tag <- basename(run_dir)
  thr_num <- suppressWarnings(as.numeric(sub("^thr_", "", thr_tag)))

  dataset_summary_path <- file.path(run_dir, "transfer", "eval", "checks", "01_dataset_summary.tsv")
  anchor_summary_path  <- file.path(run_dir, "transfer", "eval", "checks", "02_anchor_summary.tsv")
  early_summary_path   <- file.path(run_dir, "transfer", "eval", "checks", "03_early_mouse_summary.tsv")
  overall_path         <- file.path(run_dir, "transfer", "eval", "checks", "06_overall_pass_check.tsv")
  core_path            <- file.path(run_dir, "anchors", "core_gene_set.tsv")
  weights_path         <- file.path(run_dir, "anchors", "gene_weights.tsv")

  if (!file.exists(dataset_summary_path)) stop("Missing: ", dataset_summary_path)
  if (!file.exists(anchor_summary_path))  stop("Missing: ", anchor_summary_path)
  if (!file.exists(early_summary_path))   stop("Missing: ", early_summary_path)
  if (!file.exists(overall_path))         stop("Missing: ", overall_path)

  dataset_summary <- read_tsv(dataset_summary_path, show_col_types = FALSE)
  anchor_summary  <- read_tsv(anchor_summary_path, show_col_types = FALSE)
  early_summary   <- read_tsv(early_summary_path, show_col_types = FALSE)
  overall_tbl     <- read_tsv(overall_path, show_col_types = FALSE)

  n_core <- if (file.exists(core_path)) nrow(read_tsv(core_path, show_col_types = FALSE)) else NA_integer_
  n_weights <- if (file.exists(weights_path)) nrow(read_tsv(weights_path, show_col_types = FALSE)) else NA_integer_

  one_row <- tibble(
    threshold = thr_num,
    n_core_genes = n_core,
    n_weight_genes = n_weights,

    anchor_mean_delta = safe_mean(anchor_summary$mean_delta_pass),
    anchor_median_delta = safe_median(anchor_summary$mean_delta_pass),
    anchor_mean_auc = safe_mean(anchor_summary$mean_auc_pass),
    anchor_median_auc = safe_median(anchor_summary$mean_auc_pass),

    early_mean_delta = safe_mean(early_summary$mean_delta_pass),
    early_median_delta = safe_median(early_summary$mean_delta_pass),
    early_mean_auc = safe_mean(early_summary$mean_auc_pass),
    early_median_auc = safe_median(early_summary$mean_auc_pass),

    n_anchor_minimum_pass = overall_tbl$n_anchor_minimum_pass[1],
    all_true_anchors_pass_minimum = overall_tbl$all_true_anchors_pass_minimum[1],
    n_early_mouse_minimum_pass = overall_tbl$n_early_mouse_minimum_pass[1],
    prop_early_mouse_minimum_pass = overall_tbl$prop_early_mouse_minimum_pass[1]
  )

  by_dataset <- dataset_summary %>%
    mutate(threshold = thr_num) %>%
    select(
      threshold, gse_id,
      n_total_contrasts, n_pass, n_fail, pass_rate,
      mean_delta_pass, mean_auc_pass,
      prop_positive_delta_pass, prop_auc_gt_0_5_pass
    )

  list(one_row = one_row, by_dataset = by_dataset)
}

runs <- lapply(threshold_dirs, read_one_run)

summary_tbl <- bind_rows(lapply(runs, `[[`, "one_row")) %>%
  arrange(threshold)

by_dataset_tbl <- bind_rows(lapply(runs, `[[`, "by_dataset")) %>%
  arrange(gse_id, threshold)

out_summary <- file.path(sens_root, "threshold_sensitivity_summary.tsv")
out_by_dataset <- file.path(sens_root, "threshold_sensitivity_by_dataset.tsv")

write_tsv(summary_tbl, out_summary)
write_tsv(by_dataset_tbl, out_by_dataset)

message("\n==============================")
message("THRESHOLD SENSITIVITY SUMMARY")
message("==============================")
print(summary_tbl, n = nrow(summary_tbl))

message("\nWrote:")
message("  ", out_summary)
message("  ", out_by_dataset)