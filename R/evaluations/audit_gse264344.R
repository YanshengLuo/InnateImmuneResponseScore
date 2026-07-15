#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# audit_gse264344_stratified_imrs.R
#
# Purpose:
#   Stratified audit of GSE264344 IMRS behavior by:
#     1) tissue
#     2) timepoint
#     3) tissue x timepoint
#
# Inputs:
#   D:/IMRS_Project/05_score/transfer/scores/GSE264344__featurecounts_validation__imrs_scores.tsv
#   D:/IMRS_Project/00_metadata/verified_metadata/scoring/GSE264344/GSE264344_design.tsv
#
# Outputs:
#   D:/IMRS_Project/05_score/literature_compare/GSE264344_audit/
#     - GSE264344_score_plus_design.tsv
#     - GSE264344_overall_summary.tsv
#     - GSE264344_by_tissue_summary.tsv
#     - GSE264344_by_time_summary.tsv
#     - GSE264344_by_tissue_time_summary.tsv
#     - GSE264344_boxplot_tissue_time.png
#     - GSE264344_delta_by_tissue_time.png
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

score_file  <- file.path(project_root, "05_score", "transfer", "scores",
                         "GSE264344__featurecounts_validation__imrs_scores.tsv")
design_file <- file.path(project_root, "00_metadata", "verified_metadata", "scoring",
                         "GSE264344", "GSE264344_design.tsv")

out_dir <- file.path(project_root, "05_score", "literature_compare", "GSE264344_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_out      <- file.path(out_dir, "GSE264344_score_plus_design.tsv")
overall_out     <- file.path(out_dir, "GSE264344_overall_summary.tsv")
tissue_out      <- file.path(out_dir, "GSE264344_by_tissue_summary.tsv")
time_out        <- file.path(out_dir, "GSE264344_by_time_summary.tsv")
tissue_time_out <- file.path(out_dir, "GSE264344_by_tissue_time_summary.tsv")
plot1_out       <- file.path(out_dir, "GSE264344_boxplot_tissue_time.png")
plot2_out       <- file.path(out_dir, "GSE264344_delta_by_tissue_time.png")

if (!file.exists(score_file)) stop("Missing score file: ", score_file)
if (!file.exists(design_file)) stop("Missing design file: ", design_file)

# -------------------------
# HELPERS
# -------------------------
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

summarize_by_group <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols)), condition_simple) %>%
    summarise(
      n = n(),
      mean_imrs_z = safe_mean(imrs_z),
      sd_imrs_z = safe_sd(imrs_z),
      mean_imrs_raw = safe_mean(imrs_raw),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = condition_simple,
      values_from = c(n, mean_imrs_z, sd_imrs_z, mean_imrs_raw),
      names_sep = "__"
    ) %>%
    mutate(
      delta_imrs = mean_imrs_z__DELIVERY - mean_imrs_z__CONTROL
    )
}

# -------------------------
# READ AND MERGE
# -------------------------
score_df <- read_tsv(score_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    condition_simple = toupper(trimws(as.character(condition_simple))),
    imrs_raw = as.numeric(imrs_raw),
    imrs_z = as.numeric(imrs_z),
    coverage = as.numeric(coverage)
  )

design_df <- read_tsv(design_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    condition_simple = toupper(trimws(as.character(condition_simple))),
    tissue = as.character(tissue),
    timepoint_hr = as.numeric(timepoint_hr),
    time_h = as.numeric(time_h),
    group = as.character(group)
  )

merged <- score_df %>%
  left_join(
    design_df %>%
      select(sample_id, condition_simple, is_control, group, tissue, timepoint_hr, time_h, batch),
    by = c("sample_id", "condition_simple")
  ) %>%
  mutate(
    tissue = ifelse(is.na(tissue), "Unknown", tissue),
    timepoint_hr = ifelse(is.na(timepoint_hr), time_h, timepoint_hr),
    tissue_time = paste0(tissue, " | ", timepoint_hr, "h")
  )

write_tsv(merged, merged_out, na = "")

# -------------------------
# OVERALL SUMMARY
# -------------------------
overall_summary <- merged %>%
  group_by(condition_simple) %>%
  summarise(
    n = n(),
    mean_imrs_z = safe_mean(imrs_z),
    sd_imrs_z = safe_sd(imrs_z),
    mean_imrs_raw = safe_mean(imrs_raw),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = condition_simple,
    values_from = c(n, mean_imrs_z, sd_imrs_z, mean_imrs_raw),
    names_sep = "__"
  ) %>%
  mutate(
    delta_imrs = mean_imrs_z__DELIVERY - mean_imrs_z__CONTROL
  )

write_tsv(overall_summary, overall_out, na = "")

# -------------------------
# STRATIFIED SUMMARIES
# -------------------------
by_tissue <- summarize_by_group(merged, "tissue") %>%
  arrange(desc(delta_imrs), tissue)

by_time <- summarize_by_group(merged, "timepoint_hr") %>%
  arrange(timepoint_hr)

by_tissue_time <- summarize_by_group(merged, c("tissue", "timepoint_hr")) %>%
  arrange(tissue, timepoint_hr)

write_tsv(by_tissue, tissue_out, na = "")
write_tsv(by_time, time_out, na = "")
write_tsv(by_tissue_time, tissue_time_out, na = "")

# -------------------------
# PLOT 1: SAMPLE DISTRIBUTIONS
# -------------------------
plot_df <- merged %>%
  mutate(
    tissue_time = paste0(tissue, "\n", timepoint_hr, "h"),
    tissue_time = fct_inorder(tissue_time)
  )

p1 <- ggplot(plot_df, aes(x = condition_simple, y = imrs_z, color = condition_simple)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.85) +
  facet_wrap(~ tissue_time, scales = "free_y") +
  scale_color_manual(values = c("CONTROL" = "#4D4D4D", "DELIVERY" = "#B22222")) +
  labs(
    title = "GSE264344 IMRS by tissue and timepoint",
    x = "Condition",
    y = "IMRS z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(plot1_out, p1, width = 12, height = 8, dpi = 400, bg = "white")

# -------------------------
# PLOT 2: DELTA IMRS BY STRATUM
# -------------------------
plot_delta <- by_tissue_time %>%
  mutate(
    tissue = as.character(tissue),
    timepoint_hr = as.character(timepoint_hr),
    tissue_time = paste0(tissue, "\n", timepoint_hr, "h")
  ) %>%
  arrange(delta_imrs) %>%
  mutate(
    tissue_time = factor(tissue_time, levels = unique(tissue_time))
  )

p2 <- ggplot(plot_delta, aes(x = tissue_time, y = delta_imrs)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "GSE264344 delta IMRS by tissue and timepoint",
    x = "Stratum",
    y = "Delta IMRS (DELIVERY - CONTROL)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(plot2_out, p2, width = 9, height = 6, dpi = 400, bg = "white")

cat("Saved:\n")
cat("  ", merged_out, "\n")
cat("  ", overall_out, "\n")
cat("  ", tissue_out, "\n")
cat("  ", time_out, "\n")
cat("  ", tissue_time_out, "\n")
cat("  ", plot1_out, "\n")
cat("  ", plot2_out, "\n")
cat("Done.\n")