#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ============================================================
# compare_imrs_to_literature.R
#
# Purpose:
#   1. Build dataset-level IMRS summary from per-dataset score files
#   2. Compare observed IMRS behavior against paper-based expectation
#
# Inputs:
#   <project_root>/05_score/transfer/scores/*__imrs_scores.tsv
#   <project_root>/audit/imrs_literature_expectation_table_v2.tsv
#
# Optional:
#   <project_root>/audit/imrs_prototype_truth_table.tsv
#
# Outputs:
#   <project_root>/05_score/literature_compare/all_dataset_imrs_summary.tsv
#   <project_root>/05_score/literature_compare/imrs_vs_literature.tsv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

scores_dir <- file.path(project_root, "05_score", "transfer", "scores")
audit_dir  <- file.path(project_root, "audit")
out_dir    <- file.path(project_root, "05_score", "literature_compare")

lit_file   <- file.path(audit_dir, "imrs_literature_expectation_table_v2.tsv")
truth_file <- file.path(audit_dir, "imrs_prototype_truth_table.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_out <- file.path(out_dir, "all_dataset_imrs_summary.tsv")
compare_out <- file.path(out_dir, "imrs_vs_literature.tsv")

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

extract_dataset_id <- function(path) {
  b <- basename(path)
  sub("__.*$", "", b)
}

extract_score_id <- function(path) {
  sub("\\.tsv$", "", basename(path))
}

classify_observed_imrs <- function(delta_imrs, delivery_sd_imrs_z) {
  if (!is.finite(delta_imrs)) return(NA_character_)

  if (delta_imrs >= 10) {
    return("Strong")
  }

  if (is.finite(delivery_sd_imrs_z) && delivery_sd_imrs_z >= 3) {
    return("Heterogeneous")
  }

  "Weak"
}

normalize_expected_regime <- function(x) {
  x <- tolower(trimws(as.character(x)))

  case_when(
    str_detect(x, "heterogeneous") ~ "Heterogeneous",
    str_detect(x, "low") ~ "Weak",
    str_detect(x, "strong") ~ "Strong",
    TRUE ~ NA_character_
  )
}

compare_imrs_to_paper <- function(observed_regime, expected_imrs_call) {
  expected_regime <- normalize_expected_regime(expected_imrs_call)
  expected_text <- tolower(trimws(as.character(expected_imrs_call)))

  case_when(
    is.na(observed_regime) | is.na(expected_regime) ~ NA_character_,
    observed_regime == expected_regime ~ "match",
    observed_regime == "Heterogeneous" & expected_regime == "Weak" ~ "heterogeneous_vs_expected_low",
    observed_regime == "Weak" & expected_regime == "Heterogeneous" ~ "weaker_than_paper_context",
    observed_regime == "Weak" & str_detect(expected_text, "not expected") ~ "unexpectedly_low",
    observed_regime == "Heterogeneous" & str_detect(expected_text, "not expected") ~ "heterogeneous_but_paper_suggests_higher",
    observed_regime == "Strong" & expected_regime == "Weak" ~ "higher_than_expected",
    observed_regime == "Strong" & expected_regime == "Heterogeneous" ~ "stronger_than_expected",
    TRUE ~ paste0("mismatch:", observed_regime, "_vs_", expected_regime)
  )
}

read_score_file <- function(f) {
  df <- read_tsv(f, show_col_types = FALSE, progress = FALSE)

  req <- c("sample_id", "condition_simple", "imrs_z")
  missing_req <- setdiff(req, names(df))
  if (length(missing_req) > 0) {
    stop("Missing required columns in score file ", f, ": ",
         paste(missing_req, collapse = ", "))
  }

  df %>%
    mutate(
      sample_id = as.character(sample_id),
      condition_simple = toupper(trimws(as.character(condition_simple))),
      imrs_z = as.numeric(imrs_z)
    )
}

# -------------------------
# READ SCORE FILES
# -------------------------
score_files <- list.files(
  scores_dir,
  pattern = "__imrs_scores\\.tsv$",
  full.names = TRUE
)

if (length(score_files) == 0) {
  stop("No score files found in: ", scores_dir)
}

dataset_summaries <- map_dfr(score_files, function(f) {
  df <- read_score_file(f)

  ctrl <- df %>% filter(condition_simple == "CONTROL")
  del  <- df %>% filter(condition_simple == "DELIVERY")

  tibble(
    dataset_id = extract_dataset_id(f),
    id = extract_score_id(f),
    score_file = f,
    n_control = nrow(ctrl),
    n_delivery = nrow(del),
    control_mean_imrs_z = safe_mean(ctrl$imrs_z),
    control_sd_imrs_z = safe_sd(ctrl$imrs_z),
    delivery_mean_imrs_z = safe_mean(del$imrs_z),
    delivery_sd_imrs_z = safe_sd(del$imrs_z),
    delta_imrs = safe_mean(del$imrs_z) - safe_mean(ctrl$imrs_z)
  )
}) %>%
  mutate(
    observed_imrs_regime = map2_chr(
      delta_imrs, delivery_sd_imrs_z, classify_observed_imrs
    ),
    low_imrs_flag = ifelse(is.finite(delta_imrs) & delta_imrs < 10, TRUE, FALSE)
  ) %>%
  arrange(desc(delta_imrs), dataset_id)

write_tsv(dataset_summaries, summary_out, na = "")

# -------------------------
# READ LITERATURE TABLE
# -------------------------
if (!file.exists(lit_file)) {
  stop("Missing literature expectation table: ", lit_file)
}

lit_tbl <- read_tsv(lit_file, show_col_types = FALSE, progress = FALSE)

required_lit <- c(
  "dataset_id",
  "expected_imrs_call",
  "literature_comparison_take"
)
missing_lit <- setdiff(required_lit, names(lit_tbl))
if (length(missing_lit) > 0) {
  stop("Literature table must contain: ",
       paste(required_lit, collapse = ", "),
       ". Missing: ",
       paste(missing_lit, collapse = ", "))
}

compare_tbl <- dataset_summaries %>%
  left_join(lit_tbl, by = "dataset_id") %>%
  mutate(
    expected_regime = normalize_expected_regime(expected_imrs_call),
    literature_match_status = compare_imrs_to_paper(
      observed_imrs_regime, expected_imrs_call
    )
  )

# -------------------------
# OPTIONAL TRUTH TABLE
# -------------------------
if (file.exists(truth_file)) {
  truth_tbl <- read_tsv(truth_file, show_col_types = FALSE, progress = FALSE)

  label_col <- NULL
  candidates <- c("manual_label", "truth_label", "manual_class", "label")
  hit <- candidates[candidates %in% names(truth_tbl)]
  if (length(hit) > 0) label_col <- hit[1]

  if (!is.null(label_col)) {
    truth_tbl2 <- truth_tbl %>%
      transmute(
        dataset_id = as.character(dataset_id),
        truth_table_label = as.character(.data[[label_col]])
      )

    compare_tbl <- compare_tbl %>%
      left_join(truth_tbl2, by = "dataset_id")
  }
}

compare_tbl <- compare_tbl %>%
  mutate(
    audit_note = case_when(
      literature_match_status == "match" ~
        "IMRS is consistent with literature expectation.",
      literature_match_status == "unexpectedly_low" ~
        "Paper suggests stronger innate activation than observed IMRS.",
      literature_match_status == "heterogeneous_vs_expected_low" ~
        "IMRS suggests mixed or variable response rather than uniformly low biology.",
      literature_match_status == "weaker_than_paper_context" ~
        "IMRS is weaker than paper context would suggest.",
      literature_match_status == "heterogeneous_but_paper_suggests_higher" ~
        "IMRS shows heterogeneity where paper context suggests stronger response.",
      literature_match_status == "higher_than_expected" ~
        "IMRS is stronger than literature expectation.",
      literature_match_status == "stronger_than_expected" ~
        "IMRS is more strongly activated than heterogeneous paper expectation.",
      is.na(literature_match_status) ~
        "No literature expectation available.",
      TRUE ~
        "Review manually."
    )
  )

write_tsv(compare_tbl, compare_out, na = "")

cat("Saved:\n")
cat("  ", summary_out, "\n")
cat("  ", compare_out, "\n")
cat("Done.\n")