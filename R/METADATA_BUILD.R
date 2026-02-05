## ================================
## 00_metadata construction (robust)
## Dataset: <dataset_id>
## Reads: 00_metadata/<dataset_id>/SraRunTable.csv
## Writes:
##   00_metadata/<dataset_id>/<dataset_id>_samples.tsv
##   00_metadata/<dataset_id>/<dataset_id>_design.tsv   (includes tissue/timepoint/batch)
## ================================

library(readr)
library(dplyr)
library(stringr)
library(tibble)

## ----------------------------
## 1) Dataset + paths
## ----------------------------
# Allow wrapper script to set dataset_id
if (!exists("dataset_id", inherits = FALSE)) {
  stop("dataset_id is not set. Run via wrapper or set dataset_id before sourcing METADATA_BUILD.R.")
}

project_root <- "C:/Users/john/Desktop/IMRS_Project"
meta_dir <- file.path(project_root, "00_metadata", dataset_id)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

runinfo_path <- file.path(meta_dir, "SraRunTable.csv")

## ----------------------------
## 2) Helpers
## ----------------------------
pick_col <- function(df, candidates) {
  nms <- names(df)
  dn <- tolower(nms)
  cc <- tolower(candidates)
  hit <- intersect(cc, dn)
  if (length(hit) == 0) return(NA_character_)
  nms[match(hit[1], dn)]
}

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
    stop("Stopping: required metadata fields were not detected.")
  }

  cat("===============================================================\n\n")
}

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
## 4) Column synonyms
## ----------------------------
syn <- list(
  srr_id        = c("Run", "run", "SRR", "srr", "srr_id", "srr_accession"),
  biosample_id  = c("BioSample", "biosample", "BioSample_ID", "biosample_id"),
  experiment_id = c("Experiment", "experiment", "SRX", "srx", "experiment_id"),
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

required_detect <- c("srr_id", "biosample_id", "experiment_id", "timepoint", "delivery_type")
report_and_stop_if_missing(col_map, required_detect, context = paste(dataset_id, "RunInfo mapping"))

## ----------------------------
## 5) Build standardized samples table
## ----------------------------
samples <- runinfo %>%
  transmute(
    dataset_id = dataset_id,
    srr_id = as.character(.data[[col_map$srr_id]]),
    biosample_id = as.character(.data[[col_map$biosample_id]]),
    experiment_id = as.character(.data[[col_map$experiment_id]]),

    species = "Mus_musculus",

    tissue = if (!is.na(col_map$tissue)) as.character(.data[[col_map$tissue]]) else NA_character_,
    source_name = if (!is.na(col_map$source_name)) as.character(.data[[col_map$source_name]]) else NA_character_,

    timepoint_raw = as.character(.data[[col_map$timepoint]]),
    delivery_type = as.character(.data[[col_map$delivery_type]]),

    condition = NA_character_,
    control_type = NA_character_,
    payload = NA_character_,

    replicate = NA_integer_,
    batch = NA_character_,   # keep for future use / manual fill
    notes = NA_character_
  )

assert_nonempty(samples, "srr_id")
assert_nonempty(samples, "timepoint_raw")
assert_nonempty(samples, "delivery_type")

## ----------------------------
## 6) Parse timepoints → hours (NO coercion warnings)
## ----------------------------
samples <- samples %>%
  mutate(
    # normalize hidden spaces and punctuation
    timepoint_raw = str_trim(str_replace_all(timepoint_raw, "\u00A0", " ")),
    tp_txt = str_to_lower(str_replace_all(timepoint_raw, "[^a-z0-9\\.]+", " ")),
    timepoint_hr = case_when(
      tp_txt %in% c("baseline", "base line", "pre", "pretreatment", "pre treatment", "naive") ~ 0,
      str_detect(tp_txt, "\\b(h|hr|hrs|hour|hours)\\b") ~ readr::parse_number(tp_txt),
      str_detect(tp_txt, "^[0-9]+\\.?[0-9]*$") ~ readr::parse_number(tp_txt),
      TRUE ~ NA_real_
    )
  ) %>%
  select(-tp_txt)

cat("\nTimepoint diagnostics (raw unique):\n")
print(sort(unique(samples$timepoint_raw)))

na_after_parse <- samples %>%
  filter(is.na(timepoint_hr)) %>%
  count(timepoint_raw, sort = TRUE)

if (nrow(na_after_parse) > 0) {
  cat("\nWARNING: timepoint_hr is NA after parsing:\n")
  print(na_after_parse)
}

## ----------------------------
## 7) Define condition and group (dataset-specific rule)
## ----------------------------
# For GSE264344, delivery_type appears to be: baseline / Ad5 / Ad26
# Rule:
#   baseline -> control
#   everything else -> delivery
samples <- samples %>%
  mutate(
    delivery_type = str_trim(delivery_type),
    delivery_lc = str_to_lower(delivery_type),
    condition = if_else(delivery_lc == "baseline", "control", "delivery"),
    control_type = if_else(condition == "control", "naive", NA_character_),
    group = case_when(
      condition == "control"  ~ "baseline_0",
      condition == "delivery" ~ paste0("delivery_", str_replace_all(delivery_lc, "[^a-z0-9]+", "_")),
      TRUE ~ "needs_review"
    )
  ) %>%
  select(-delivery_lc)

cat("\nCondition assignment counts:\n")
print(table(samples$condition, useNA = "ifany"))

cat("\nGroup counts:\n")
print(table(samples$group, useNA = "ifany"))

## ----------------------------
## 8) Final checks
## ----------------------------
cat("\nBaseline controls by tissue/timepoint:\n")
print(samples %>% filter(condition == "control") %>% count(tissue, timepoint_hr, sort = TRUE))

cat("\nDelivery samples by tissue/delivery/timepoint:\n")
print(samples %>% filter(condition == "delivery") %>% count(tissue, delivery_type, timepoint_hr, sort = TRUE))

if (!any(samples$group == "baseline_0")) {
  stop("No baseline_0 samples detected. DE contrasts against baseline_0 will fail.")
}

if (any(is.na(samples$timepoint_hr))) {
  stop("Some samples have timepoint_hr = NA. Fix parsing rules before proceeding.")
}

## ----------------------------
## 9) Save outputs
## ----------------------------
samples_path <- file.path(meta_dir, paste0(dataset_id, "_samples.tsv"))
write_tsv(samples, samples_path)

# design.tsv now includes tissue, timepoint_hr, batch as requested
design_out <- samples %>%
  transmute(
    sample_id = srr_id,     # must match count matrix column names
    group = group,
    tissue = tissue,
    timepoint_hr = timepoint_hr,
    batch = batch
  )

design_path <- file.path(meta_dir, paste0(dataset_id, "_design.tsv"))
write_tsv(design_out, design_path)

cat("\nWrote:\n")
cat("  -", samples_path, "\n")
cat("  -", design_path, "\n\n")
