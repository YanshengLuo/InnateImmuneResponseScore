suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# =========================
# Config
# =========================
imrs_root <- "D:/IMRS_Project"
framework_root <- "D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/Dataset_Level_Confidence_score_prototype"

design_dir <- file.path(imrs_root, "00_metadata", "verified_metadata", "scoring")
rawmeta_dir <- file.path(imrs_root, "00_metadata", "raw_metadata")

out_file <- file.path(framework_root, "inputs", "02_sample_metadata.tsv")
log_file <- file.path(framework_root, "logs", "02_build_sample_metadata.log")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

# =========================
# Logging helper
# =========================
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

# =========================
# Helpers
# =========================
extract_dataset_id_from_name <- function(x) {
  stringr::str_extract(x, "GSE\\d+")
}

safe_read_design <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_rawmeta_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_rawmeta_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

# =========================
# Find files
# =========================
design_files <- list.files(
  design_dir,
  pattern = "_design\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

rawmeta_files <- character(0)
if (dir.exists(rawmeta_dir)) {
  rawmeta_files <- c(
    list.files(rawmeta_dir, pattern = "SraRunTable.*\\.csv$", recursive = TRUE, full.names = TRUE),
    list.files(rawmeta_dir, pattern = "SraRunTable.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  )
}

log_message("Design dir: ", design_dir)
log_message("Raw metadata dir: ", rawmeta_dir)
log_message("Found design files: ", length(design_files))
log_message("Found raw metadata files: ", length(rawmeta_files))

if (length(design_files) > 0) {
  log_message("Example design files:")
  for (x in head(design_files, 5)) log_message("  ", x)
}

if (length(rawmeta_files) > 0) {
  log_message("Example raw metadata files:")
  for (x in head(rawmeta_files, 5)) log_message("  ", x)
}

if (length(design_files) == 0) {
  stop("No design files found in: ", design_dir)
}

# =========================
# Build raw metadata lookup by Run
# =========================
rawmeta_lookup_list <- purrr::map(rawmeta_files, function(f) {
  df <- tryCatch(
    if (grepl("\\.csv$", f, ignore.case = TRUE)) safe_read_rawmeta_csv(f) else safe_read_rawmeta_tsv(f),
    error = function(e) NULL
  )
  if (is.null(df)) {
    log_message("Failed to read raw metadata: ", f)
    return(NULL)
  }

  run_col <- names(df)[tolower(names(df)) == "run"]
  if (length(run_col) != 1) {
    log_message("No Run column found in raw metadata file: ", f)
    return(NULL)
  }

  get_col <- function(nm) {
    hit <- names(df)[tolower(names(df)) == tolower(nm)]
    if (length(hit) == 1) return(df[[hit]])
    rep(NA, nrow(df))
  }

  tibble(
    sample_id = as.character(df[[run_col]]),
    platform = as.character(get_col("Platform")),
    library_layout = as.character(get_col("LibraryLayout")),
    library_source = as.character(get_col("LibrarySource")),
    raw_treatment = as.character(get_col("treatment"))
  )
})

rawmeta_lookup <- dplyr::bind_rows(rawmeta_lookup_list)

if (nrow(rawmeta_lookup) == 0) {
  rawmeta_lookup <- tibble(
    sample_id = character(),
    platform = character(),
    library_layout = character(),
    library_source = character(),
    raw_treatment = character()
  )
} else {
  rawmeta_lookup <- rawmeta_lookup %>%
    filter(!is.na(sample_id), sample_id != "") %>%
    distinct(sample_id, .keep_all = TRUE)
}

log_message("Raw metadata sample lookup rows: ", nrow(rawmeta_lookup))

# =========================
# Read verified design files
# =========================
sample_metadata_list <- purrr::map(design_files, function(f) {
  dataset_id <- extract_dataset_id_from_name(f)
  if (is.na(dataset_id)) {
    log_message("Skipping design file with no dataset ID: ", f)
    return(NULL)
  }

  df <- safe_read_design(f)

  if (!("sample_id" %in% names(df))) stop("Missing sample_id in: ", f)
  if (!("group" %in% names(df))) stop("Missing group in: ", f)

  if (!("tissue" %in% names(df))) df$tissue <- NA_character_
  if (!("timepoint_hr" %in% names(df))) df$timepoint_hr <- NA_real_
  if (!("batch" %in% names(df))) df$batch <- NA_character_

  df %>%
    transmute(
      dataset_id = dataset_id,
      sample_id = as.character(sample_id),
      group = as.character(group),
      tissue = as.character(tissue),
      time_hr = suppressWarnings(as.numeric(timepoint_hr)),
      batch = as.character(batch)
    ) %>%
    group_by(dataset_id, group) %>%
    mutate(replicate = row_number()) %>%
    ungroup()
})

sample_metadata <- dplyr::bind_rows(sample_metadata_list)

if (nrow(sample_metadata) == 0) {
  stop("No sample metadata could be built from design files.")
}

# =========================
# Optional enrichment from raw metadata
# =========================
sample_metadata <- sample_metadata %>%
  left_join(rawmeta_lookup, by = "sample_id") %>%
  arrange(dataset_id, group, replicate, sample_id)

write_tsv(sample_metadata, out_file)
log_message("Saved sample metadata: ", out_file)
log_message("Rows written: ", nrow(sample_metadata))