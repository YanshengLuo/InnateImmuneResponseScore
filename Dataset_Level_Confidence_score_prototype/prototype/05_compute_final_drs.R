#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# ============================================================
# 05_compute_final_drs.R
#
# Purpose:
#   Aggregate calibrated component scores into final prototype DRS
#
# Inputs:
#   intermediate/04_drs_component_metrics.tsv
#
# Outputs:
#   outputs/05_final_drs.tsv
#   outputs/05_final_drs_ranked.tsv
#   outputs/05_final_drs_summary.tsv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
framework_root <- if (length(args) >= 1) args[1] else
  "D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/Dataset_Level_Confidence_score_prototype"

intermediate_dir <- file.path(framework_root, "intermediate")
outputs_dir <- file.path(framework_root, "outputs")
logs_dir <- file.path(framework_root, "logs")

dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

component_file <- file.path(intermediate_dir, "04_drs_component_metrics.tsv")
out_file <- file.path(outputs_dir, "05_final_drs.tsv")
ranked_file <- file.path(outputs_dir, "05_final_drs_ranked.tsv")
summary_file <- file.path(outputs_dir, "05_final_drs_summary.tsv")
log_file <- file.path(logs_dir, "05_compute_final_drs.log")

if (file.exists(log_file)) file.remove(log_file)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

assign_drs_class <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 0.80 ~ "High",
    x >= 0.60 ~ "Moderate",
    x >= 0.40 ~ "Low",
    TRUE ~ "Very Low"
  )
}

stopifnot(file.exists(component_file))

df <- read_tsv(component_file, show_col_types = FALSE, progress = FALSE)

required_cols <- c("dataset_id", "Ccoh", "Cout", "Csep", "Cctx")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in 04_drs_component_metrics.tsv: ",
       paste(missing_cols, collapse = ", "))
}

final_drs <- df %>%
  mutate(
    DRS = (Ccoh + Cout + Csep + Cctx) / 4,
    DRS_class = assign_drs_class(DRS)
  ) %>%
  select(
    dataset_id,
    dataset_class,
    organism,
    assay_type,
    tissue,
    treatment_type,
    n_samples_used,
    n_genes_selected,
    Ccoh,
    Cout,
    Csep,
    Cctx,
    primary_failure_mode,
    DRS,
    DRS_class,
    everything()
  )

ranked_drs <- final_drs %>%
  arrange(desc(DRS)) %>%
  mutate(DRS_rank = row_number()) %>%
  select(DRS_rank, everything())

summary_tbl <- ranked_drs %>%
  summarise(
    n_datasets = n(),
    mean_DRS = mean(DRS, na.rm = TRUE),
    median_DRS = median(DRS, na.rm = TRUE),
    min_DRS = min(DRS, na.rm = TRUE),
    max_DRS = max(DRS, na.rm = TRUE)
  )

write_tsv(final_drs, out_file, na = "")
write_tsv(ranked_drs, ranked_file, na = "")
write_tsv(summary_tbl, summary_file, na = "")

log_message("Saved final DRS: ", out_file)
log_message("Saved ranked DRS: ", ranked_file)
log_message("Saved DRS summary: ", summary_file)
log_message("Datasets aggregated: ", nrow(final_drs))
log_message("Done.")