#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
framework_root <- if (length(args) >= 1) args[1] else
  "D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/Dataset_Level_Confidence_score_prototype"

# CHANGE THIS if your IMRS summary file lives somewhere else
imrs_summary_file <- if (length(args) >= 2) args[2] else
  "D:/IMRS_Project/05_score/failure_diagnosis/all_dataset_imrs_summary.tsv"

drs_rank_file <- file.path(framework_root, "outputs", "05_final_drs_ranked.tsv")
out_file <- file.path(framework_root, "outputs", "06_drs_vs_imrs.tsv")

stopifnot(file.exists(drs_rank_file))
stopifnot(file.exists(imrs_summary_file))

drs <- read_tsv(drs_rank_file, show_col_types = FALSE, progress = FALSE)
imrs <- read_tsv(imrs_summary_file, show_col_types = FALSE, progress = FALSE)

# adjust these if needed depending on your IMRS table columns
required_imrs <- c("dataset_id")
missing_imrs <- setdiff(required_imrs, names(imrs))
if (length(missing_imrs) > 0) {
  stop("IMRS file missing required columns: ", paste(missing_imrs, collapse = ", "))
}

merged <- drs %>%
  left_join(imrs, by = "dataset_id")

write_tsv(merged, out_file, na = "")
cat("Saved:", out_file, "\n")