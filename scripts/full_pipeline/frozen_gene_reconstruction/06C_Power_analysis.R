#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 06C — Power / information screening (ANCHOR ONLY)
# - Uses se_meta from Step 6B output: gene_heterogeneity.tsv
# - Computes power to detect |log2FC| = Delta at alpha=0.05 (two-sided)
# - Computes for ALL genes, but flags low-power only for core genes (if core exists)
#
# Inputs:
#   <project_root>/05_score/anchors/gene_heterogeneity.tsv
#   <project_root>/05_score/anchors/core_gene_set.tsv   (optional; used for flagging)
#
# Outputs:
#   <project_root>/05_score/anchors/gene_power.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# ---- Phase fixed to anchor only (consistent with your plan) ----
PHASE <- "anchor"

phase_out_dir <- function(phase) {
  if (phase == "anchor") return("anchors")
  if (phase == "calibration") return("calibration")
  stop("bad phase: ", phase)
}

# Best-practice defaults for your use case
Delta <- 1.0         # detectable effect size in log2FC units (absolute)
alpha <- 0.05        # two-sided
power_cutoff <- 0.80 # standard

# Numerical safety
se_floor <- 1e-8

# Two-sided normal-approx power for |Z| > z_alpha/2 under true mean = Delta/SE
power_two_sided <- function(se, Delta = 1.0, alpha = 0.05) {
  se <- pmax(se, se_floor)
  zcrit <- qnorm(1 - alpha/2)
  delta_z <- Delta / se
  # Power = P(Z > zcrit) + P(Z < -zcrit), for Z ~ N(delta_z, 1)
  p_upper <- 1 - pnorm(zcrit - delta_z)
  p_lower <- pnorm(-zcrit - delta_z)
  p_upper + p_lower
}

# -------------------------
# MAIN
# -------------------------
out_root <- file.path(project_root, "05_score", phase_out_dir(PHASE))
if (!dir.exists(out_root)) {
  stop("Missing folder: ", out_root, "\nDid Step 06A/06B write outputs here?")
}

het_path <- file.path(out_root, "gene_heterogeneity.tsv")
if (!file.exists(het_path)) {
  stop("Missing Step 06B output: ", het_path)
}

het <- read_tsv(het_path, show_col_types = FALSE)

if (!all(c("gene_id", "se_meta") %in% names(het))) {
  stop("gene_heterogeneity.tsv must contain columns: gene_id, se_meta\nFile: ", het_path)
}

# Core genes (optional) used only for flagging
core_path <- file.path(out_root, "core_gene_set.tsv")
core_genes <- character(0)
if (file.exists(core_path)) {
  core_df <- read_tsv(core_path, show_col_types = FALSE)
  if ("gene_id" %in% names(core_df)) core_genes <- as.character(core_df$gene_id)
}

out <- het %>%
  transmute(
    gene_id   = as.character(.data$gene_id),
    phase     = PHASE,
    se_meta   = suppressWarnings(as.numeric(.data$se_meta)),
    beta_meta = if ("beta_meta" %in% names(het)) suppressWarnings(as.numeric(.data$beta_meta)) else NA_real_
  ) %>%
  mutate(
    is_core = .data$gene_id %in% core_genes,
    power_delta = ifelse(is.na(.data$se_meta), NA_real_,
                         power_two_sided(.data$se_meta, Delta = Delta, alpha = alpha)),
    low_power_flag = .data$is_core & !is.na(.data$power_delta) & (.data$power_delta < power_cutoff)
  ) %>%
  arrange(desc(.data$is_core), desc(.data$low_power_flag), .data$power_delta, .data$gene_id)

out_path <- file.path(out_root, "gene_power.tsv")
write_tsv(out, out_path)

message("PHASE=", PHASE, ": wrote ", out_path,
        " | genes=", nrow(out),
        " | core_loaded=", length(core_genes),
        " | low_power_core=", sum(out$low_power_flag, na.rm = TRUE),
        " | Delta=", Delta, " | alpha=", alpha, " | cutoff=", power_cutoff)

message("\nDONE: Step 06C power computed using SE_meta (two-sided normal approx).")

