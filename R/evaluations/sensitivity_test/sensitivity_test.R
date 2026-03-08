#!/usr/bin/env Rscript
# ============================================================
# IMRS threshold sensitivity driver
#
# Runs:
#   06A -> 06B -> 06C -> 07 -> 08 -> 09 -> 09_check
#
# for multiple effect-size thresholds and snapshots outputs.
#
# REQUIREMENT:
#   06A_core_gene_set.R must accept:
#     Rscript 06A_core_gene_set.R <project_root> <effect_threshold>
#
# Assumes other scripts already work with:
#     Rscript <script> <project_root>
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"
scripts_root <- if (length(args) >= 2) args[2] else file.path(project_root, "Hypergator_scripts", "InnateImmuneResponseScore", "R")

thresholds <- c(0.4, 0.5, 0.6)

# ---- script paths ----
script_06A <- file.path(scripts_root, "06A_core_gene_set.R")
script_06B <- file.path(scripts_root, "06B_gene_heterogeneity.R")
script_06C <- file.path(scripts_root, "06C_Power_analysis.R")
script_07  <- file.path(scripts_root, "07_gene_weights.R")
script_08  <- file.path(scripts_root, "08_score_samples.R")
script_09  <- file.path(scripts_root, "09_split_eval.R")
script_09check <- file.path(scripts_root, "09_check_eval_summary.R")

needed_scripts <- c(script_06A, script_06B, script_06C, script_07, script_08, script_09, script_09check)
missing_scripts <- needed_scripts[!file.exists(needed_scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing scripts:\n", paste(missing_scripts, collapse = "\n"))
}

# ---- output base for snapshots ----
sens_root <- file.path(project_root, "05_score", "sensitivity_runs")
dir.create(sens_root, recursive = TRUE, showWarnings = FALSE)

# ---- helper: run system command and stop on failure ----
run_rscript <- function(script, args_vec = character()) {
  cmd <- c(shQuote(script), vapply(args_vec, shQuote, character(1)))
  full_cmd <- paste("Rscript", paste(cmd, collapse = " "))
  message("\n[RUN] ", full_cmd)
  status <- system(full_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ":\n", full_cmd)
  }
}

# ---- helper: safe copy file if exists ----
copy_if_exists <- function(from, to) {
  if (file.exists(from)) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE)
    TRUE
  } else {
    FALSE
  }
}

# ---- helper: safe copy directory recursively ----
copy_dir_if_exists <- function(from, to) {
  if (dir.exists(from)) {
    if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, recursive = TRUE)
    TRUE
  } else {
    FALSE
  }
}

# ---- helper: snapshot outputs after each threshold run ----
snapshot_outputs <- function(thr) {
  thr_tag <- paste0("thr_", format(thr, nsmall = 1))
  out_dir <- file.path(sens_root, thr_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Anchor model outputs
  anchor_dir <- file.path(project_root, "05_score", "anchors")
  copy_dir_if_exists(anchor_dir, file.path(out_dir, "anchors"))

  # Mouse transfer/eval outputs
  transfer_dir <- file.path(project_root, "05_score", "transfer")
  copy_dir_if_exists(transfer_dir, file.path(out_dir, "transfer"))

  # Human transfer outputs (optional; usually unaffected by threshold only if weights changed, so keep snapshot)
  human_dir <- file.path(project_root, "05_score", "human_transfer")
  copy_dir_if_exists(human_dir, file.path(out_dir, "human_transfer"))

  # Figures from Step09 mouse eval
  fig_dir <- file.path(project_root, "06_figures", "step09_calibration")
  copy_dir_if_exists(fig_dir, file.path(out_dir, "figures", "step09_calibration"))

  message("[SNAPSHOT] Saved outputs to: ", out_dir)
}

# -------------------------
# MAIN LOOP
# -------------------------
for (thr in thresholds) {
  message("\n====================================================")
  message("THRESHOLD SENSITIVITY RUN: ", thr)
  message("====================================================")

  # Step 06A: threshold-sensitive
  run_rscript(script_06A, c(project_root, as.character(thr)))

  # Step 06B: uses anchor outputs from 06A
  run_rscript(script_06B, c(project_root))

  # Step 06C: uses heterogeneity output from 06B
  run_rscript(script_06C, c(project_root))

  # Step 07: rebuild weights from current anchor outputs
  run_rscript(script_07, c(project_root))

  # Step 08: rescore mouse datasets with new frozen weights
  run_rscript(script_08, c(project_root))

  # Step 09: reevaluate mouse contrasts
  run_rscript(script_09, c(project_root))

  # Step 09 post-check summary
  run_rscript(script_09check, c(project_root))

  # Snapshot outputs from this threshold run
  snapshot_outputs(thr)
}

message("\nDONE: threshold sensitivity runs completed.")
message("Snapshots written under: ", sens_root)