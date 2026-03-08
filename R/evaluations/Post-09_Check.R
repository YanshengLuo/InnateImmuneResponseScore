#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 09 CHECKS — Post-evaluation summary and pass/fail
#
# Reads:
#   <project_root>/05_score/transfer/eval/step09_split_eval.tsv
#
# Writes:
#   <project_root>/05_score/transfer/eval/checks/
#       01_dataset_summary.tsv
#       02_anchor_summary.tsv
#       03_early_mouse_summary.tsv
#       04_anchor_pass_check.tsv
#       05_early_mouse_pass_check.tsv
#       06_overall_pass_check.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

if (!requireNamespace("pROC", quietly = TRUE)) {
  message("pROC not found. Installing...")
  install.packages("pROC", repos = "https://cloud.r-project.org")
}

library(pROC)
# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# TRUE training anchors only
anchor_ids <- c(
  "GSE39129",
  "GSE167521",
  "GSE264344"
)

# All mouse datasets in the <=24 h evaluation pool
early_mouse_ids <- c(
  "GSE39129",
  "GSE167521",
  "GSE264344",
  "GSE279372",
  "GSE279744",
  "GSE262515"
)

eval_path <- file.path(project_root, "05_score", "transfer", "eval", "step09_split_eval.tsv")
out_dir   <- file.path(project_root, "05_score", "transfer", "eval", "checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(eval_path)) {
  stop("Missing Step 09 eval table: ", eval_path)
}

# -------------------------
# HELPERS
# -------------------------
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

safe_prop <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

# -------------------------
# READ EVAL TABLE
# -------------------------
eval_tbl <- read_tsv(eval_path, show_col_types = FALSE)

required_cols <- c("gse_id", "pass", "delta_mean_imrs_z", "auc_imrs_z")
miss <- setdiff(required_cols, names(eval_tbl))
if (length(miss) > 0) {
  stop("Eval table missing required columns: ", paste(miss, collapse = ", "))
}

# -------------------------
# 1) DATASET-LEVEL SUMMARY
# -------------------------
dataset_summary <- eval_tbl %>%
  group_by(gse_id) %>%
  summarise(
    n_total_contrasts = n(),
    n_pass = sum(pass, na.rm = TRUE),
    n_fail = sum(!pass, na.rm = TRUE),
    pass_rate = n_pass / n_total_contrasts,

    mean_delta_all = safe_mean(delta_mean_imrs_z),
    median_delta_all = safe_median(delta_mean_imrs_z),
    mean_auc_all = safe_mean(auc_imrs_z),

    mean_delta_pass = safe_mean(delta_mean_imrs_z[pass]),
    median_delta_pass = safe_median(delta_mean_imrs_z[pass]),
    mean_auc_pass = safe_mean(auc_imrs_z[pass]),

    prop_positive_delta_pass = safe_prop(delta_mean_imrs_z[pass] > 0),
    prop_auc_gt_0_5_pass = safe_prop(auc_imrs_z[pass] > 0.5),

    .groups = "drop"
  ) %>%
  arrange(desc(n_pass), desc(mean_delta_pass))

# -------------------------
# 2) TRUE ANCHOR SUMMARY
# -------------------------
anchor_summary <- eval_tbl %>%
  filter(gse_id %in% anchor_ids) %>%
  group_by(gse_id) %>%
  summarise(
    n_total_contrasts = n(),
    n_pass = sum(pass, na.rm = TRUE),
    n_fail = sum(!pass, na.rm = TRUE),
    pass_rate = n_pass / n_total_contrasts,

    mean_delta_pass = safe_mean(delta_mean_imrs_z[pass]),
    median_delta_pass = safe_median(delta_mean_imrs_z[pass]),
    mean_auc_pass = safe_mean(auc_imrs_z[pass]),

    prop_positive_delta_pass = safe_prop(delta_mean_imrs_z[pass] > 0),
    prop_auc_gt_0_5_pass = safe_prop(auc_imrs_z[pass] > 0.5),

    .groups = "drop"
  ) %>%
  arrange(match(gse_id, anchor_ids))

anchor_pass_check <- anchor_summary %>%
  mutate(
    has_passing_contrasts = n_pass > 0,
    mean_delta_positive = is.finite(mean_delta_pass) & mean_delta_pass > 0,
    mean_auc_above_random = is.finite(mean_auc_pass) & mean_auc_pass > 0.5,
    anchor_pass_minimum = has_passing_contrasts & mean_delta_positive & mean_auc_above_random
  )

# -------------------------
# 3) ALL EARLY MOUSE SUMMARY
# -------------------------
early_mouse_summary <- eval_tbl %>%
  filter(gse_id %in% early_mouse_ids) %>%
  group_by(gse_id) %>%
  summarise(
    n_total_contrasts = n(),
    n_pass = sum(pass, na.rm = TRUE),
    n_fail = sum(!pass, na.rm = TRUE),
    pass_rate = n_pass / n_total_contrasts,

    mean_delta_pass = safe_mean(delta_mean_imrs_z[pass]),
    median_delta_pass = safe_median(delta_mean_imrs_z[pass]),
    mean_auc_pass = safe_mean(auc_imrs_z[pass]),

    prop_positive_delta_pass = safe_prop(delta_mean_imrs_z[pass] > 0),
    prop_auc_gt_0_5_pass = safe_prop(auc_imrs_z[pass] > 0.5),

    .groups = "drop"
  ) %>%
  arrange(match(gse_id, early_mouse_ids))

early_mouse_pass_check <- early_mouse_summary %>%
  mutate(
    has_passing_contrasts = n_pass > 0,
    mean_delta_positive = is.finite(mean_delta_pass) & mean_delta_pass > 0,
    mean_auc_above_random = is.finite(mean_auc_pass) & mean_auc_pass > 0.5,
    early_mouse_pass_minimum = has_passing_contrasts & mean_delta_positive & mean_auc_above_random
  )

# -------------------------
# 4) OVERALL PASS CHECK
# -------------------------
overall_pass_check <- tibble(
  n_anchor_datasets = length(anchor_ids),
  n_anchor_present = sum(anchor_ids %in% unique(eval_tbl$gse_id)),
  n_anchor_minimum_pass = sum(anchor_pass_check$anchor_pass_minimum, na.rm = TRUE),
  all_true_anchors_pass_minimum = all(anchor_pass_check$anchor_pass_minimum, na.rm = TRUE),

  n_early_mouse_datasets = length(early_mouse_ids),
  n_early_mouse_present = sum(early_mouse_ids %in% unique(eval_tbl$gse_id)),
  n_early_mouse_minimum_pass = sum(early_mouse_pass_check$early_mouse_pass_minimum, na.rm = TRUE),
  prop_early_mouse_minimum_pass = mean(early_mouse_pass_check$early_mouse_pass_minimum, na.rm = TRUE)
)

# -------------------------
# WRITE TABLES
# -------------------------
write_tsv(dataset_summary, file.path(out_dir, "01_dataset_summary.tsv"))
write_tsv(anchor_summary, file.path(out_dir, "02_anchor_summary.tsv"))
write_tsv(early_mouse_summary, file.path(out_dir, "03_early_mouse_summary.tsv"))
write_tsv(anchor_pass_check, file.path(out_dir, "04_anchor_pass_check.tsv"))
write_tsv(early_mouse_pass_check, file.path(out_dir, "05_early_mouse_pass_check.tsv"))
write_tsv(overall_pass_check, file.path(out_dir, "06_overall_pass_check.tsv"))

# -------------------------
# PRINT
# -------------------------
message("\n==============================")
message("DATASET SUMMARY")
message("==============================")
print(dataset_summary, n = nrow(dataset_summary))

message("\n==============================")
message("TRUE ANCHOR SUMMARY")
message("==============================")
print(anchor_summary, n = nrow(anchor_summary))

message("\n==============================")
message("ANCHOR PASS CHECK")
message("==============================")
print(anchor_pass_check, n = nrow(anchor_pass_check))

message("\n==============================")
message("EARLY MOUSE SUMMARY (<=24 h pool)")
message("==============================")
print(early_mouse_summary, n = nrow(early_mouse_summary))

message("\n==============================")
message("EARLY MOUSE PASS CHECK")
message("==============================")
print(early_mouse_pass_check, n = nrow(early_mouse_pass_check))

message("\n==============================")
message("OVERALL PASS CHECK")
message("==============================")
print(overall_pass_check)

message("\nWrote checks to: ", out_dir)