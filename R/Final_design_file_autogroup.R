# --------------------------------------------
# Use code below to group your design.tsv after verifing the accruacy of your data manually.
# In this code, we separated all samples as replicates with same tisse/timepoint/delivery_method/dosage
# Build per-stratum (tissue x timepoint) design.tsv files
# and keep tissue/timepoint columns (redundant by request)
#
# Inputs:
#   1) meta_tsv: original metadata TSV (your "original .tsv")
#   2) out_dir: output folder for new design files
#
# Assumptions:
#   - meta has columns that can be mapped to: sample_id, tissue, timepoint_hr, delivery_type
#   - controls are identified by delivery_type matching control_keywords
#   - batch is optional; if present and you set use_batch=TRUE, it will stratify on batch too
# --------------------------------------------

# ==== User settings ====
meta_tsv <- "metadata/GSEXXXXX/metadata.tsv"   # <-- change
gse_id   <- "GSEXXXXX"                         # <-- change
out_dir  <- file.path("metadata", gse_id, "design_strata")

# How to recognize controls (edit as needed)
control_keywords <- c("control", "pbs", "vehicle", "saline", "naive", "untreated", "mock")

# Batch behavior
use_batch <- TRUE  # set FALSE if you want to ignore batch even if present

# Minimum replicates (per group) to keep a contrast
min_reps_per_group <- 2

# ==== Libraries ====
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# ==== Helpers ====
normalize_str <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower()
}

is_control <- function(delivery, control_keywords) {
  normalize_str(delivery) %in% control_keywords
}

sanitize_token <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_sub(1, 80)
}

# Column mapping: tries to find likely columns if your names differ
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ==== Read metadata ====
meta <- read_tsv(meta_tsv, show_col_types = FALSE)

# Try to auto-map common column names (edit if needed)
col_sample <- pick_col(meta, c("sample_id","sample","Sample","SampleID","run","Run","SRR","srr_id","srr_accession"))
col_tissue <- pick_col(meta, c("tissue","Tissue","organ","Organ","tissue_type"))
col_time   <- pick_col(meta, c("timepoint_hr","timepoint","time_point","Time","time","hours_post","hrs_post"))
col_deliv  <- pick_col(meta, c("delivery","delivery_type","platform","vector","treatment","group","arm","condition"))

col_batch  <- pick_col(meta, c("batch","Batch","lane","Lane","flowcell","Flowcell","run_batch","seq_batch"))

stopifnot(!is.na(col_sample), !is.na(col_tissue), !is.na(col_time), !is.na(col_deliv))

# Standardize columns we will use
df <- meta %>%
  transmute(
    sample_id   = .data[[col_sample]],
    tissue      = .data[[col_tissue]],
    timepoint_hr = .data[[col_time]],
    delivery_raw = .data[[col_deliv]],
    batch       = if (!is.na(col_batch)) .data[[col_batch]] else NA_character_
  ) %>%
  mutate(
    tissue       = as.character(tissue),
    timepoint_hr = as.character(timepoint_hr),
    delivery_raw = as.character(delivery_raw),
    batch        = as.character(batch),
    delivery_norm = normalize_str(delivery_raw),
    is_control    = is_control(delivery_norm, control_keywords),
    delivery      = ifelse(is_control, "control", delivery_raw) # keep original delivery labels for treated
  )

# Create output folder
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define stratum keys
if (use_batch && !all(is.na(df$batch))) {
  df <- df %>% mutate(stratum_batch = ifelse(is.na(batch) | batch == "", "NA", batch))
  strata <- df %>% distinct(tissue, timepoint_hr, stratum_batch)
} else {
  df <- df %>% mutate(stratum_batch = "NA")
  strata <- df %>% distinct(tissue, timepoint_hr, stratum_batch)
}

# ==== Main loop ====
written_files <- list()
dropped_strata <- list()

for (i in seq_len(nrow(strata))) {
  tissue_i <- strata$tissue[i]
  time_i   <- strata$timepoint_hr[i]
  batch_i  <- strata$stratum_batch[i]

  sub <- df %>%
    filter(
      tissue == tissue_i,
      timepoint_hr == time_i,
      stratum_batch == batch_i
    )

  # Must have at least one control and one treated group
  if (sum(sub$is_control, na.rm = TRUE) == 0 || sum(!sub$is_control, na.rm = TRUE) == 0) {
    dropped_strata[[length(dropped_strata) + 1]] <- tibble(
      tissue = tissue_i, timepoint_hr = time_i, batch = batch_i,
      reason = "Missing control or missing treated samples in this stratum"
    )
    next
  }

  # Enforce minimum replicates per group (control + each treated group)
  counts <- sub %>%
    count(delivery, name = "n")

  # control reps
  n_control <- counts %>% filter(delivery == "control") %>% pull(n)
  if (length(n_control) == 0) n_control <- 0

  # if control < min
  if (n_control < min_reps_per_group) {
    dropped_strata[[length(dropped_strata) + 1]] <- tibble(
      tissue = tissue_i, timepoint_hr = time_i, batch = batch_i,
      reason = paste0("Control replicates < ", min_reps_per_group)
    )
    next
  }

  # Identify treated groups meeting min reps; keep only those + control
  ok_treated <- counts %>%
    filter(delivery != "control", n >= min_reps_per_group) %>%
    pull(delivery)

  if (length(ok_treated) == 0) {
    dropped_strata[[length(dropped_strata) + 1]] <- tibble(
      tissue = tissue_i, timepoint_hr = time_i, batch = batch_i,
      reason = paste0("No treated group has >= ", min_reps_per_group, " replicates")
    )
    next
  }

  sub2 <- sub %>%
    filter(delivery == "control" | delivery %in% ok_treated) %>%
    mutate(
      # redundant columns requested
      tissue = tissue_i,
      timepoint_hr = time_i
    ) %>%
    select(sample_id, tissue, timepoint_hr, batch, delivery)

  # Output filename: design_GSE_tissue_timeh(_batchX).tsv
  tissue_tag <- sanitize_token(tissue_i)
  time_tag   <- sanitize_token(paste0(time_i, "h"))
  batch_tag  <- if (!identical(batch_i, "NA")) paste0("_batch_", sanitize_token(batch_i)) else ""

  out_file <- file.path(out_dir, paste0("design_", gse_id, "_", tissue_tag, "_", time_tag, batch_tag, ".tsv"))

  write_tsv(sub2, out_file)
  written_files[[length(written_files) + 1]] <- tibble(
    file = out_file, tissue = tissue_i, timepoint_hr = time_i, batch = batch_i,
    n_samples = nrow(sub2),
    n_control = sum(sub2$delivery == "control"),
    treated_groups = paste(unique(sub2$delivery[sub2$delivery != "control"]), collapse = ";")
  )
}

# ==== Summary outputs ====
written_df <- if (length(written_files)) bind_rows(written_files) else tibble()
dropped_df <- if (length(dropped_strata)) bind_rows(dropped_strata) else tibble()

message("Wrote ", nrow(written_df), " design files to: ", out_dir)
if (nrow(dropped_df) > 0) {
  message("Dropped ", nrow(dropped_df), " strata (see dropped_strata.tsv).")
}

# Save summaries
if (nrow(written_df) > 0) write_tsv(written_df, file.path(out_dir, "written_files.tsv"))
if (nrow(dropped_df) > 0) write_tsv(dropped_df, file.path(out_dir, "dropped_strata.tsv"))

# Print a quick view
written_df
