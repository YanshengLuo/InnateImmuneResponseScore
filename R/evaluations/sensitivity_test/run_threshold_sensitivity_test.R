#!/usr/bin/env Rscript
# ============================================================
# IMRS threshold sensitivity driver
# Corrected for your actual file names / directories
# with robust directory snapshot copying via fs::dir_copy()
# ============================================================

suppressPackageStartupMessages({
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

modeling_root <- file.path(project_root, "Hypergator_scripts", "InnateImmuneResponseScore", "R", "modeling")
eval_root     <- file.path(project_root, "Hypergator_scripts", "InnateImmuneResponseScore", "R", "evaluations")

thresholds <- c(0.4, 0.5, 0.6)

# -------------------------
# actual script paths
# -------------------------
script_06A <- file.path(modeling_root, "06A_core_gene_set.R")
script_06B <- file.path(modeling_root, "06B_gene_heterogeneity.R")
script_06C <- file.path(modeling_root, "06C_Power_analysis.R")
script_07  <- file.path(modeling_root, "07_weight_estimation.R")
script_08  <- file.path(modeling_root, "08_score_samples.R")
script_09  <- file.path(modeling_root, "09_calibration_evaluation.R")
script_09check <- file.path(eval_root, "Post-09_Check.R")

needed_scripts <- c(
  script_06A, script_06B, script_06C,
  script_07, script_08, script_09, script_09check
)

missing_scripts <- needed_scripts[!file.exists(needed_scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing scripts:\n", paste(missing_scripts, collapse = "\n"))
}

# -------------------------
# output snapshot root
# -------------------------
sens_root <- file.path(project_root, "05_score", "sensitivity_runs")
dir.create(sens_root, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# helpers
# -------------------------
run_rscript <- function(script, args_vec = character()) {
  cmd <- c(shQuote(script), vapply(args_vec, shQuote, character(1)))
  full_cmd <- paste("Rscript", paste(cmd, collapse = " "))
  message("\n[RUN] ", full_cmd)
  status <- system(full_cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ":\n", full_cmd)
  }
}

copy_dir_if_exists <- function(from, to) {
  if (!dir.exists(from)) {
    warning("Source directory does not exist: ", from)
    return(FALSE)
  }

  if (dir.exists(to)) {
    unlink(to, recursive = TRUE, force = TRUE)
  }

  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)

  fs::dir_copy(from, to, overwrite = TRUE)
  TRUE
}

snapshot_outputs <- function(thr) {
  thr_tag <- paste0("thr_", format(thr, nsmall = 1))
  out_dir <- file.path(sens_root, thr_tag)

  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  copy_dir_if_exists(
    file.path(project_root, "05_score", "anchors"),
    file.path(out_dir, "anchors")
  )

  copy_dir_if_exists(
    file.path(project_root, "05_score", "transfer"),
    file.path(out_dir, "transfer")
  )

  if (dir.exists(file.path(project_root, "05_score", "human_transfer"))) {
    copy_dir_if_exists(
      file.path(project_root, "05_score", "human_transfer"),
      file.path(out_dir, "human_transfer")
    )
  }

  if (dir.exists(file.path(project_root, "06_figures", "step09_calibration"))) {
    copy_dir_if_exists(
      file.path(project_root, "06_figures", "step09_calibration"),
      file.path(out_dir, "figures", "step09_calibration")
    )
  }

  message("[SNAPSHOT] Saved outputs to: ", out_dir)
}

# -------------------------
# MAIN LOOP
# -------------------------
for (thr in thresholds) {
  message("\n====================================================")
  message("THRESHOLD SENSITIVITY RUN: ", thr)
  message("====================================================")

  # IMPORTANT:
  # 06A must accept <project_root> <effect_threshold>
  run_rscript(script_06A, c(project_root, as.character(thr)))

  run_rscript(script_06B, c(project_root))
  run_rscript(script_06C, c(project_root))
  run_rscript(script_07,  c(project_root))
  run_rscript(script_08,  c(project_root))
  run_rscript(script_09,  c(project_root))
  run_rscript(script_09check, c(project_root))

  snapshot_outputs(thr)
}

message("\nDONE: threshold sensitivity runs completed.")
message("Snapshots written under: ", sens_root)