#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 00X — Build DATASET-LEVEL design.tsv files for SCORING (Step 08)
#
# INPUT (your verified metadata folder):
#   <project_root>/00_metadata/verified_metadata/*.tsv
#   (e.g., GSE39129_design.tsv, GSE167521_design.tsv, ...)
#
# OUTPUT (scoring subfolder, dataset-level):
#   <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#
# Also writes:
#   <project_root>/00_metadata/verified_metadata/scoring/control_audit.tsv
#
# Notes:
# - This uses the SAME control-detection logic pattern you used for splitting:
#     baseline_0 -> known_controls_by_file -> keyword fallback
# - Adds: is_control + condition_simple (CONTROL/DELIVERY)
# - Dataset-level design includes all samples across tissues/time/batch.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ----------------------------
# USER SETTINGS (Windows)
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# IMPORTANT: your corrected directory
in_dir <- file.path(project_root, "00_metadata", "verified_metadata")

# Output into the same verified_metadata folder, under scoring/
out_root <- file.path(in_dir, "scoring")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", project_root)
message("Input dir:     ", in_dir)
message("Output dir:    ", out_root)

# ----------------------------
# Helpers (mirrors split script style)
# ----------------------------
norm_label <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

safe_token <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- str_trim(x)
  x <- str_replace_all(x, "[\\\\/:*?\"<>|]+", "_")
  x <- str_replace_all(x, "\\s+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_+|_+$", "")
  ifelse(nchar(x) == 0, "NA", x)
}

to_hours <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "h$", "")
  suppressWarnings(as.numeric(x))
}

# ----------------------------
# Control detection (same logic as your split script)
# ----------------------------
known_controls_by_file <- list(
  "GSE190850_HUMAN_design.tsv" = c("delivery_unstimulated"),
  "GSE39129_design.tsv"        = c("delivery_none")
)

control_keywords <- c(
  "baseline", "control", "ctrl",
  "pbs", "saline", "vehicle",
  "mock", "sham", "empty",
  "naive", "untreated", "uninfected",
  "unstimulated", "none", "no_treatment", "pre"
)
control_regex <- paste0("(", paste(unique(norm_label(control_keywords)), collapse = "|"), ")")

detect_controls <- function(df, file_base) {
  groups <- unique(as.character(df$group))
  controls <- character(0)

  # Priority 1: baseline_0 if present
  if ("baseline_0" %in% groups) controls <- c(controls, "baseline_0")

  # Priority 2: file-specific known control labels
  if (file_base %in% names(known_controls_by_file)) {
    controls <- c(controls, known_controls_by_file[[file_base]])
  }

  controls <- unique(controls)
  if (length(controls) > 0) return(controls)

  # Priority 3: keyword fallback
  gnorm <- norm_label(groups)
  groups[str_detect(gnorm, control_regex)] %>% unique()
}

# ----------------------------
# Discover TSVs in verified_metadata
# ----------------------------
if (!dir.exists(in_dir)) stop("Input dir does not exist: ", in_dir)

tsv_files <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)

# Exclude derived artifacts and split outputs (same as your splitter exclusions + scoring dir)
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__split_index\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__contrast_index\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__group_map\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "^control_report\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "split_design")]
tsv_files <- tsv_files[!str_detect(dirname(tsv_files), paste0("[/\\\\]", "splited", "([/\\\\]|$)"))]
tsv_files <- tsv_files[!str_detect(dirname(tsv_files), paste0("[/\\\\]", "scoring", "([/\\\\]|$)"))]

if (length(tsv_files) == 0) stop("No eligible TSV files found in: ", in_dir)

message("Found ", length(tsv_files), " input TSV files.")

# ----------------------------
# MAIN
# ----------------------------
audit_rows <- list()
written_n <- 0L
skipped_n <- 0L

for (f in tsv_files) {
  file_base <- basename(f)
  stem <- tools::file_path_sans_ext(file_base) # e.g., GSE39129_design
  message("\n=== Processing: ", file_base, " ===")

  df <- read_tsv(f, show_col_types = FALSE)

  req <- c("sample_id", "group", "tissue", "timepoint_hr", "batch")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    message("SKIP (missing columns): ", paste(miss, collapse = ", "))
    skipped_n <- skipped_n + 1L
    next
  }

  # Infer dataset ID from filename prefix before "_design"
  dataset_id <- str_match(stem, "^(.+?)_design")[, 2]
  if (is.na(dataset_id) || dataset_id == "") {
    dataset_id <- str_match(stem, "^(GSE[^_]+)")[, 2]
  }
  if (is.na(dataset_id) || dataset_id == "") {
    message("SKIP (cannot infer dataset_id from filename): ", file_base)
    skipped_n <- skipped_n + 1L
    next
  }

  df2 <- df %>%
    mutate(
      sample_id    = trimws(as.character(sample_id)),
      group        = trimws(as.character(group)),
      tissue       = as.character(tissue),
      timepoint_hr = as.character(timepoint_hr),
      time_h       = to_hours(timepoint_hr),
      batch        = ifelse(is.na(batch) | batch == "", "NA", as.character(batch))
    )

  if (any(!nzchar(df2$sample_id))) stop("Empty sample_id detected in: ", file_base)

  control_labels <- detect_controls(df2, file_base)

  if (length(control_labels) == 0) {
    message("WARN: No controls detected in ", file_base, ". Skipping dataset scoring design.")
    audit_rows[[length(audit_rows) + 1]] <- tibble(
      dataset_id = dataset_id,
      file = file_base,
      n_samples = nrow(df2),
      n_controls = 0L,
      controls = "",
      status = "SKIP_no_controls"
    )
    skipped_n <- skipped_n + 1L
    next
  } else {
    message("Controls detected: ", paste(control_labels, collapse = ", "))
  }

  df2 <- df2 %>%
    mutate(
      is_control = group %in% control_labels,
      condition_simple = ifelse(is_control, "CONTROL", "DELIVERY")
    )

  n_ctrl <- sum(df2$is_control)
  n_del  <- sum(!df2$is_control)

  # Deduplicate sample_id if needed (keep first; warn)
  if (any(duplicated(df2$sample_id))) {
    dup_ids <- unique(df2$sample_id[duplicated(df2$sample_id)])
    message("WARN: duplicated sample_id(s) in ", dataset_id, ": ",
            paste(head(dup_ids, 10), collapse = ", "),
            if (length(dup_ids) > 10) " ..." else "")
    df2 <- df2 %>% distinct(sample_id, .keep_all = TRUE)
  }

  # Write dataset-level design into verified_metadata/scoring/<DATASET>/
  ds_dir <- file.path(out_root, safe_token(dataset_id))
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)

  out_path <- file.path(ds_dir, paste0(safe_token(dataset_id), "_design.tsv"))

  # Stable column order (Step 08 minimum + useful fields)
  keep_first <- c("sample_id", "condition_simple", "is_control",
                  "group", "tissue", "timepoint_hr", "time_h", "batch")
  keep_first <- keep_first[keep_first %in% names(df2)]

  df_out <- df2 %>% select(all_of(keep_first), everything())

  write_tsv(df_out, out_path)
  message("WROTE scoring design: ", out_path)
  written_n <- written_n + 1L

  audit_rows[[length(audit_rows) + 1]] <- tibble(
    dataset_id = dataset_id,
    file = file_base,
    out_path = out_path,
    n_samples = nrow(df_out),
    n_controls = as.integer(n_ctrl),
    n_delivery = as.integer(n_del),
    controls = paste(control_labels, collapse = ";"),
    status = ifelse(n_ctrl > 0 && n_del > 0, "OK", "WARN_only_one_class")
  )
}

audit <- bind_rows(audit_rows)
audit_path <- file.path(out_root, "control_audit.tsv")
write_tsv(audit, audit_path)

message("\nDONE.")
message("Wrote scoring designs: ", written_n)
message("Skipped: ", skipped_n)
message("Audit: ", audit_path)