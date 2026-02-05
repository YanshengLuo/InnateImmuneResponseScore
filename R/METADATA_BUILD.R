## ================================
## 00_metadata construction 
## Dataset: <dataset_id>
## Reads: 00_metadata/<dataset_id>/SraRunTable.csv
## Writes: 00_metadata/<dataset_id>/<dataset_id>_samples.tsv
##         00_metadata/<dataset_id>/<dataset_id>_design.tsv
## ================================

library(readr)
library(dplyr)
library(stringr)
library(tibble)

## ----------------------------
## 1) Dataset + paths
## ----------------------------
dataset_id <- "GSE264344"

project_root <- "C:/Users/john/Desktop/IMRS_Project"
meta_dir <- file.path(project_root, "00_metadata", dataset_id)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

runinfo_path <- file.path(meta_dir, "SraRunTable.csv")

## ----------------------------
## 2) Helpers: column detection + validation reporting
## ----------------------------

# pick the first matching column (case-insensitive)
pick_col <- function(df, candidates) {
  nms <- names(df)
  dn <- tolower(nms)
  cc <- tolower(candidates)
  hit <- intersect(cc, dn)
  if (length(hit) == 0) return(NA_character_)
  nms[match(hit[1], dn)]
}

# print a report and stop if required columns are missing
report_and_stop_if_missing <- function(col_map, required, context = "") {
  cat("\n================ METADATA COLUMN DETECTION REPORT ================\n")
  if (nzchar(context)) cat("Context:", context, "\n\n")

  for (nm in names(col_map)) {
    cat(sprintf("  %-14s -> %s\n", nm, ifelse(is.na(col_map[[nm]]), "NOT FOUND", col_map[[nm]])))
  }

  missing_required <- required[is.na(unlist(col_map[required]))]

  cat("\nRequired fields:\n")
  if (length(missing_required) == 0) {
    cat("  OK: all required fields detected.\n")
  } else {
    cat("  MISSING required:", paste(missing_required, collapse = ", "), "\n")
    cat("===============================================================\n\n")
    stop("Stopping: required metadata fields were not detected. Update synonyms or source file.")
  }

  cat("===============================================================\n\n")
}

# require non-empty values in critical fields
assert_nonempty <- function(df, col) {
  bad <- is.na(df[[col]]) | trimws(as.character(df[[col]])) == ""
  if (any(bad)) {
    bad_n <- sum(bad)
    stop(sprintf("Column '%s' has %d missing/empty values. Fix metadata before proceeding.", col, bad_n))
  }
}

## ----------------------------
## 3) Read RunInfo
## ----------------------------
runinfo <- read_csv(runinfo_path, show_col_types = FALSE)

cat("\nLoaded RunInfo:\n")
cat("  File:", runinfo_path, "\n")
cat("  Rows:", nrow(runinfo), "  Cols:", ncol(runinfo), "\n")

## ----------------------------
## 4) Column synonym sets (edit once, reuse forever)
## ----------------------------
syn <- list(
  srr_id        = c("Run", "run", "SRR", "srr", "srr_id", "srr_accession"),
  biosample_id  = c("BioSample", "biosample", "BioSample_ID", "biosample_id"),
  experiment_id = c("Experiment", "experiment", "SRX", "srx", "experiment_id"),

  # these vary across sources
  tissue        = c("tissue", "organ", "Organ", "tissue_type"),
  source_name   = c("source_name", "source_name_ch1", "SourceName", "source"),
  timepoint     = c("timepoint", "time_point", "time", "Time", "collection_time", "hours_post", "hrs_post"),
  delivery_type = c("vaccine_type", "vaccine", "treatment", "treat", "platform", "vector",
                    "delivery", "group", "arm")
)

col_map <- list(
  srr_id        = pick_col(runinfo, syn$srr_id),
  biosample_id  = pick_col(runinfo, syn$biosample_id),
  experiment_id = pick_col(runinfo, syn$experiment_id),
  tissue        = pick_col(runinfo, syn$tissue),
  source_name   = pick_col(runinfo, syn$source_name),
  timepoint     = pick_col(runinfo, syn$timepoint),
  delivery_type = pick_col(runinfo, syn$delivery_type)
)

# Required for this dataset construction
required_detect <- c("srr_id", "biosample_id", "experiment_id", "timepoint", "delivery_type")
report_and_stop_if_missing(col_map, required_detect, context = paste(dataset_id, "RunInfo mapping"))

## ----------------------------
## 5) Initialize standardized samples table
## ----------------------------
samples <- runinfo %>%
  transmute(
    dataset_id = dataset_id,
    srr_id = as.character(.data[[col_map$srr_id]]),
    biosample_id = as.character(.data[[col_map$biosample_id]]),
    experiment_id = as.character(.data[[col_map$experiment_id]]),

    species = "Mus_musculus",   # confirmed from paper

    # optional columns: if not found, fill NA but do not crash
    tissue = if (!is.na(col_map$tissue)) as.character(.data[[col_map$tissue]]) else NA_character_,
    source_name = if (!is.na(col_map$source_name)) as.character(.data[[col_map$source_name]]) else NA_character_,

    timepoint_raw = as.character(.data[[col_map$timepoint]]),
    delivery_type = as.character(.data[[col_map$delivery_type]]),

    condition = NA_character_,      # filled later
    control_type = NA_character_,   # filled later
    payload = NA_character_,        # not needed for this dataset

    replicate = NA_integer_,
    batch = NA_character_,
    notes = NA_character_
  )

# critical field sanity
assert_nonempty(samples, "srr_id")
assert_nonempty(samples, "timepoint_raw")
assert_nonempty(samples, "delivery_type")

## ----------------------------
## 6) Parse timepoints (→ hours) [robust]
## ----------------------------
# Handles: "baseline", "0", "4h", "4 hr", "24 hours", numeric strings
samples <- samples %>%
  mutate(
    timepoint_raw = str_trim(timepoint_raw),
    timepoint_txt = str_to_lower(str_replace_all(timepoint_raw, "[^a-z0-9\\.]+", " ")),
    timepoint_hr = case_when(
      timepoint_txt %in% c("baseline", "base line", "pre", "pretreatment", "pre treatment", "naive") ~ 0,
      str_detect(timepoint_txt, "\\b[0-9]+\\.?[0-9]*\\s*h\\b|\\b[0-9]+\\.?[0-9]*\\s*hr\\b|\\b[0-9]+\\.?[0-9]*\\s*hour\\b") ~
        as.numeric(str_extract(timepoint_txt, "[0-9]+\\.?[0-9]*")),
      str_detect(timepoint_txt, "^[0-9]+\\.?[0-9]*$") ~ as.numeric(timepoint_txt),
      TRUE ~ NA_real_
    )
  ) %>%
  select(-timepoint_txt)

cat("\nTimepoint diagnostics (raw unique):\n")
print(sort(unique(samples$timepoint_raw)))

na_tp <- samples %>% filter(is.na(timepoint_hr)) %>% count(timepoint_raw, sort = TRUE)
if (nrow(na_tp) > 0) {
  cat("\nWARNING: Unparsed timepoints (timepoint_hr = NA):\n")
  print(na_tp)
}

## ----------------------------
## 7) Define condition and control (keyword rules)
## ----------------------------
baseline_kw <- c("baseline", "naive", "control", "mock", "pbs", "saline", "vehicle", "untreated", "pre")
delivery_kw <- c("vaccine", "vaccinated", "treatment", "treated", "ad26", "ad5", "adenovirus",
                 "vector", "dose", "lnp", "mrna", "sirna")

samples <- samples %>%
  mutate(
    delivery_type = str_trim(delivery_type),
    delivery_txt = str_to_lower(str_replace_all(delivery_type, "[^a-z0-9]+", " ")),
    condition = case_when(
      str_detect(delivery_txt, str_c("\\b(", str_c(baseline_kw, collapse="|"), ")\\b")) ~ "control",
      str_detect(delivery_txt, str_c("\\b(", str_c(delivery_kw, collapse="|"), ")\\b")) ~ "delivery",
      TRUE ~ "needs_review"
    ),
    control_type = if_else(condition == "control", "naive", NA_character_)
  ) %>%
  select(-delivery_txt)

cat("\nCondition assignment counts:\n")
print(table(samples$condition, useNA = "ifany"))

if (any(samples$condition == "needs_review")) {
  cat("\nWARNING: Some samples could not be confidently classified (needs_review).\n")
  cat("First 10 examples (srr_id, delivery_type):\n")
  print(samples %>% filter(condition == "needs_review") %>% select(srr_id, delivery_type) %>% head(10))
}

# Create DESeq2-ready group label
samples <- samples %>%
  mutate(
    group = case_when(
      condition == "control"  ~ "baseline_0",
      condition == "delivery" ~ paste0("delivery_", str_replace_all(str_to_lower(delivery_type), "[^a-z0-9]+", "_")),
      TRUE ~ "needs_review"
    )
  )

cat("\nGroup counts:\n")
print(table(samples$group, useNA = "ifany"))

## ----------------------------
## 8) Final structural checks
## ----------------------------

cat("\nBaseline controls by tissue/timepoint:\n")
print(
  samples %>%
    filter(condition == "control") %>%
    count(tissue, timepoint_hr, sort = TRUE)
)

cat("\nDelivery samples by tissue/delivery/timepoint:\n")
print(
  samples %>%
    filter(condition == "delivery") %>%
    count(tissue, delivery_type, timepoint_hr, sort = TRUE)
)

# Hard stop if no baseline_0 exists (DE will fail)
if (!any(samples$group == "baseline_0")) {
  stop("No baseline_0 samples detected. DE contrasts against baseline_0 will fail.")
}

# Hard stop if any needs_review remain (recommended)
if (any(samples$group == "needs_review")) {
  stop("Some samples are 'needs_review'. Fix classification rules or add overrides before proceeding.")
}

## ----------------------------
## 9) Save outputs to 00_metadata/<dataset_id>/
## ----------------------------

# full standardized metadata (rich)
samples_path <- file.path(meta_dir, paste0(dataset_id, "_samples.tsv"))
write_tsv(samples, samples_path)

# minimal DESeq2 design matrix (what your DE scripts need)
# IMPORTANT: sample_id must match count matrix column names (often SRR ids)
design_out <- samples %>%
  transmute(
    sample_id = srr_id,
    group = group
  )

design_path <- file.path(meta_dir, paste0(dataset_id, "_design.tsv"))
write_tsv(design_out, design_path)

cat("\nWrote:\n")
cat("  -", samples_path, "\n")
cat("  -", design_path, "\n\n")
