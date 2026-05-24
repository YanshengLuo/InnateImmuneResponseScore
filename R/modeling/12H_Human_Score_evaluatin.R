#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 09H — Human transfer evaluation
#
# Purpose:
#   Evaluate mouse->human transferred IMRS on the scored human dataset.
#   NO retraining, NO recalibration, NO weight updates.
#
# Inputs:
#   1) Step 08H human score file:
#      <project_root>/05_score/human_transfer/scores/*__imrs_scores.tsv
#
# Outputs:
#   Tables:
#     <project_root>/05_score/human_transfer/eval/step09_human_eval.tsv
#     <project_root>/05_score/human_transfer/eval/step09_human_sample_level.tsv
#     <project_root>/05_score/human_transfer/eval/step09_human_roc.tsv   (optional)
#
#   Figures:
#     <project_root>/06_figures/step09_human_transfer/
#
# Primary metric:
#   delta_mean_imrs_z = mean(DELIVERY) - mean(CONTROL)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# exact dataset folder/prefix, e.g. GSE190850_HUMAN
dataset_id <- if (length(args) >= 2) args[2] else "GSE190850_HUMAN"

do_auc <- TRUE
write_roc <- TRUE

scores_dir <- file.path(project_root, "05_score", "human_transfer", "scores")
eval_dir   <- file.path(project_root, "05_score", "human_transfer", "eval")
fig_dir    <- file.path(project_root, "06_figures", "step09_human_transfer")

dir.create(eval_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(scores_dir)) stop("Missing scores_dir: ", scores_dir)

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

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

cohens_d <- function(x_treat, x_ctrl) {
  x_treat <- x_treat[is.finite(x_treat)]
  x_ctrl  <- x_ctrl[is.finite(x_ctrl)]

  if (length(x_treat) < 2 || length(x_ctrl) < 2) return(NA_real_)

  m1 <- mean(x_treat)
  m0 <- mean(x_ctrl)
  s1 <- sd(x_treat)
  s0 <- sd(x_ctrl)

  sp <- sqrt(
    (((length(x_treat) - 1) * s1^2) + ((length(x_ctrl) - 1) * s0^2)) /
      (length(x_treat) + length(x_ctrl) - 2)
  )

  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m0) / sp
}

compute_auc <- function(df) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    return(list(auc = NA_real_, roc_tbl = NULL, auc_note = "pROC_not_installed"))
  }

  tmp <- df %>%
    filter(is.finite(imrs_z), !is.na(label))

  if (nrow(tmp) == 0 || length(unique(tmp$label)) < 2) {
    return(list(auc = NA_real_, roc_tbl = NULL, auc_note = "one_class_only"))
  }

  roc_obj <- pROC::roc(
    response = tmp$label,
    predictor = tmp$imrs_z,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_tbl <- tibble(
    threshold   = roc_obj$thresholds,
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )

  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    roc_tbl = roc_tbl,
    auc_note = "ok"
  )
}

save_plot_both <- function(plot_obj, out_stub, width = 8, height = 6) {
  ggsave(
    filename = paste0(out_stub, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    filename = paste0(out_stub, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
}

choose_best_score_file <- function(x) {
  if (length(x) == 1) return(x)
  preferred <- x[grepl("__counts__imrs_scores\\.tsv$", basename(x))]
  if (length(preferred) >= 1) return(preferred[1])
  x[1]
}

extract_score_id <- function(path) {
  sub("__imrs_scores\\.tsv$", "", basename(path))
}

# -------------------------
# FIND HUMAN SCORE FILE
# -------------------------
all_score_files <- list.files(
  scores_dir,
  pattern = "__imrs_scores\\.tsv$",
  full.names = TRUE
)

if (length(all_score_files) == 0) {
  stop("No Step 08H score files found under: ", scores_dir)
}

matched_files <- all_score_files[grepl(paste0("^", dataset_id, "__"), basename(all_score_files))]

if (length(matched_files) == 0) {
  stop("No human score file found for dataset_id: ", dataset_id)
}

score_path <- choose_best_score_file(matched_files)
score_id <- extract_score_id(score_path)

message("Using score file: ", score_path)

# -------------------------
# READ SCORES
# -------------------------
scores_df <- read_tsv(score_path, show_col_types = FALSE, progress = FALSE)

required_score_cols <- c("sample_id", "condition_simple", "imrs_raw", "imrs_z")
miss_score <- setdiff(required_score_cols, names(scores_df))
if (length(miss_score) > 0) {
  stop("Step 08H score file missing columns: ", paste(miss_score, collapse = ", "))
}

scores_df <- scores_df %>%
  mutate(
    sample_id = trimws(as.character(sample_id)),
    condition_simple = toupper(trimws(as.character(condition_simple))),
    imrs_raw = suppressWarnings(as.numeric(imrs_raw)),
    imrs_z   = suppressWarnings(as.numeric(imrs_z))
  ) %>%
  filter(condition_simple %in% c("CONTROL", "DELIVERY"))

if (nrow(scores_df) == 0) {
  stop("No usable rows after filtering CONTROL/DELIVERY in: ", score_path)
}

ctrl <- scores_df %>% filter(condition_simple == "CONTROL")
del  <- scores_df %>% filter(condition_simple == "DELIVERY")

n_ctrl <- nrow(ctrl)
n_del  <- nrow(del)

if (n_ctrl < 1 || n_del < 1) {
  stop("Need both CONTROL and DELIVERY samples. Found n_ctrl=", n_ctrl, ", n_del=", n_del)
}

# -------------------------
# METRICS
# -------------------------
ctrl_mean   <- safe_mean(ctrl$imrs_z)
del_mean    <- safe_mean(del$imrs_z)
ctrl_median <- safe_median(ctrl$imrs_z)
del_median  <- safe_median(del$imrs_z)

delta_mean   <- del_mean - ctrl_mean
delta_median <- del_median - ctrl_median
d_val        <- cohens_d(del$imrs_z, ctrl$imrs_z)

t_p <- tryCatch(
  t.test(del$imrs_z, ctrl$imrs_z)$p.value,
  error = function(e) NA_real_
)

wilcox_p <- tryCatch(
  wilcox.test(del$imrs_z, ctrl$imrs_z)$p.value,
  error = function(e) NA_real_
)

auc_val <- NA_real_
auc_note <- NA_character_
roc_tbl <- NULL

if (isTRUE(do_auc)) {
  auc_input <- scores_df %>%
    mutate(label = ifelse(condition_simple == "DELIVERY", 1L, 0L))

  auc_out <- compute_auc(auc_input)
  auc_val <- auc_out$auc
  auc_note <- auc_out$auc_note
  roc_tbl <- auc_out$roc_tbl
}

# -------------------------
# OUTPUT TABLES
# -------------------------
eval_tbl <- tibble(
  dataset_id = dataset_id,
  score_id = score_id,
  pass = TRUE,

  n_samples = nrow(scores_df),
  n_controls = n_ctrl,
  n_delivery = n_del,

  ortholog_coverage_fixed = if ("ortholog_coverage_fixed" %in% names(scores_df)) {
    safe_mean(scores_df$ortholog_coverage_fixed)
  } else {
    NA_real_
  },

  scoring_coverage = if ("scoring_coverage" %in% names(scores_df)) {
    safe_mean(scores_df$scoring_coverage)
  } else if ("coverage" %in% names(scores_df)) {
    safe_mean(scores_df$coverage)
  } else {
    NA_real_
  },

  control_mean_imrs_z = ctrl_mean,
  delivery_mean_imrs_z = del_mean,
  delta_mean_imrs_z = delta_mean,

  control_median_imrs_z = ctrl_median,
  delivery_median_imrs_z = del_median,
  delta_median_imrs_z = delta_median,

  control_sd_imrs_z = safe_sd(ctrl$imrs_z),
  delivery_sd_imrs_z = safe_sd(del$imrs_z),

  cohens_d_imrs_z = d_val,
  t_test_p_imrs_z = t_p,
  wilcox_p_imrs_z = wilcox_p,

  auc_imrs_z = auc_val,
  auc_note = auc_note
)

sample_tbl <- scores_df %>%
  transmute(
    dataset_id = dataset_id,
    score_id = score_id,
    sample_id = sample_id,
    condition_simple = condition_simple,
    imrs_raw = imrs_raw,
    imrs_z = imrs_z,
    ortholog_coverage_fixed = if ("ortholog_coverage_fixed" %in% names(scores_df)) ortholog_coverage_fixed else NA_real_,
    scoring_coverage = if ("scoring_coverage" %in% names(scores_df)) scoring_coverage else if ("coverage" %in% names(scores_df)) coverage else NA_real_
  )

out_eval    <- file.path(eval_dir, "step09_human_eval.tsv")
out_samples <- file.path(eval_dir, "step09_human_sample_level.tsv")
out_roc     <- file.path(eval_dir, "step09_human_roc.tsv")

write_tsv(eval_tbl, out_eval)
write_tsv(sample_tbl, out_samples)

if (isTRUE(write_roc) && !is.null(roc_tbl)) {
  write_tsv(roc_tbl, out_roc)
}

# -------------------------
# PLOTS
# -------------------------
# 1. IMRS_z by condition
p1 <- ggplot(scores_df, aes(x = condition_simple, y = imrs_z)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, height = 0) +
  theme_bw() +
  labs(
    title = paste0("Human transfer IMRS_z by condition: ", dataset_id),
    x = "Condition",
    y = "IMRS_z"
  )
save_plot_both(p1, file.path(fig_dir, "01_imrsz_by_condition"), 7, 6)

# 2. IMRS_z density
p2 <- ggplot(scores_df, aes(x = imrs_z, linetype = condition_simple)) +
  geom_density() +
  theme_bw() +
  labs(
    title = paste0("Human transfer IMRS_z density: ", dataset_id),
    x = "IMRS_z",
    y = "Density"
  )
save_plot_both(p2, file.path(fig_dir, "02_imrsz_density"), 8, 6)

# 3. IMRS_raw by condition
p3 <- ggplot(scores_df, aes(x = condition_simple, y = imrs_raw)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, height = 0) +
  theme_bw() +
  labs(
    title = paste0("Human transfer IMRS_raw by condition: ", dataset_id),
    x = "Condition",
    y = "IMRS_raw"
  )
save_plot_both(p3, file.path(fig_dir, "03_imrsraw_by_condition"), 7, 6)

# 4. ROC curve if available
if (!is.null(roc_tbl) && nrow(roc_tbl) > 0) {
  roc_plot_df <- roc_tbl %>%
    mutate(
      fpr = 1 - specificity,
      tpr = sensitivity
    )

  p4 <- ggplot(roc_plot_df, aes(x = fpr, y = tpr)) +
    geom_path() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_bw() +
    labs(
      title = paste0("Human transfer ROC: ", dataset_id, " (AUC=", round(auc_val, 3), ")"),
      x = "False positive rate",
      y = "True positive rate"
    )
  save_plot_both(p4, file.path(fig_dir, "04_roc_curve"), 7, 6)
}

# 5. Ranked samples
rank_df <- scores_df %>%
  arrange(imrs_z) %>%
  mutate(sample_rank = row_number())

p5 <- ggplot(rank_df, aes(x = sample_rank, y = imrs_z, shape = condition_simple)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = paste0("Human transfer ranked samples: ", dataset_id),
    x = "Ranked sample index",
    y = "IMRS_z"
  )
save_plot_both(p5, file.path(fig_dir, "05_ranked_samples"), 8, 6)

# -------------------------
# FINISH
# -------------------------
message("\nDONE: Step 09H human transfer evaluation finished.")
message("Eval table:   ", out_eval)
message("Sample table: ", out_samples)
if (isTRUE(write_roc) && !is.null(roc_tbl)) {
  message("ROC table:    ", out_roc)
}
message("Figures dir:  ", fig_dir)

message("\nHuman transfer result summary:")
print(eval_tbl)