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

manual_file <- file.path(framework_root, "inputs", "dataset_annotation_manual.tsv")
out_file <- file.path(framework_root, "inputs", "01_dataset_manifest.tsv")
log_file <- file.path(framework_root, "logs", "01_build_dataset_manifest.log")

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

extract_single_or_mixed <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  if (length(x) == 1) return(x)
  "mixed"
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

infer_assay_type <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  if (all(grepl("rna-seq|expression profiling by high throughput sequencing|transcript", x))) {
    return("bulk_rnaseq")
  }
  extract_single_or_mixed(x)
}

infer_organism <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x1 <- unique(x)
  if (length(x1) > 1) return("mixed")
  val <- tolower(x1)
  if (grepl("mus musculus|mouse", val)) return("mouse")
  if (grepl("homo sapiens|human", val)) return("human")
  x1
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
if (!file.exists(manual_file)) {
  stop("Missing manual annotation file: ", manual_file)
}

# =========================
# Summarize verified design files
# =========================
design_summary_list <- purrr::map(design_files, function(f) {
  dataset_id <- extract_dataset_id_from_name(f)
  if (is.na(dataset_id)) {
    log_message("Skipping design file with no dataset ID: ", f)
    return(NULL)
  }

  df <- safe_read_design(f)

  tibble(
    dataset_id = dataset_id,
    tissue_design = if ("tissue" %in% names(df)) extract_single_or_mixed(df$tissue) else NA_character_,
    time_hr_design = if ("timepoint_hr" %in% names(df)) extract_single_or_mixed(df$timepoint_hr) else NA_character_,
    grouping_variable = if ("group" %in% names(df)) "group" else NA_character_,
    n_samples_design = nrow(df)
  )
})

design_summary <- dplyr::bind_rows(design_summary_list)

if (nrow(design_summary) == 0) {
  design_summary <- tibble(
    dataset_id = character(),
    tissue_design = character(),
    time_hr_design = character(),
    grouping_variable = character(),
    n_samples_design = numeric()
  )
} else {
  design_summary <- design_summary %>%
    distinct(dataset_id, .keep_all = TRUE)
}

# =========================
# Summarize raw metadata files
# =========================
rawmeta_summary_list <- purrr::map(rawmeta_files, function(f) {
  dataset_id <- extract_dataset_id_from_name(f)
  if (is.na(dataset_id)) {
    log_message("Raw metadata file has no dataset ID in name/path, skipping manifest enrichment: ", f)
    return(NULL)
  }

  df <- tryCatch(
    if (grepl("\\.csv$", f, ignore.case = TRUE)) safe_read_rawmeta_csv(f) else safe_read_rawmeta_tsv(f),
    error = function(e) NULL
  )
  if (is.null(df)) {
    log_message("Failed to read raw metadata: ", f)
    return(NULL)
  }

  tissue_col <- names(df)[tolower(names(df)) == "tissue"]
  organism_col <- names(df)[tolower(names(df)) == "organism"]
  assay_col <- names(df)[tolower(names(df)) %in% c("assay type", "assay_type")]
  treatment_col <- names(df)[tolower(names(df)) == "treatment"]

  tibble(
    dataset_id = dataset_id,
    tissue_raw = if (length(tissue_col) == 1) extract_single_or_mixed(df[[tissue_col]]) else NA_character_,
    organism_raw = if (length(organism_col) == 1) infer_organism(df[[organism_col]]) else NA_character_,
    assay_type_raw = if (length(assay_col) == 1) infer_assay_type(df[[assay_col]]) else NA_character_,
    raw_treatment = if (length(treatment_col) == 1) extract_single_or_mixed(df[[treatment_col]]) else NA_character_
  )
})

rawmeta_summary <- dplyr::bind_rows(rawmeta_summary_list)

if (nrow(rawmeta_summary) == 0) {
  rawmeta_summary <- tibble(
    dataset_id = character(),
    tissue_raw = character(),
    organism_raw = character(),
    assay_type_raw = character(),
    raw_treatment = character()
  )
} else {
  rawmeta_summary <- rawmeta_summary %>%
    distinct(dataset_id, .keep_all = TRUE)
}

# =========================
# Manual annotation
# =========================
manual_df <- read_tsv(manual_file, show_col_types = FALSE, progress = FALSE)

required_manual <- c(
  "dataset_id", "dataset_class", "organism", "assay_type",
  "treatment_type", "include_in_analysis", "notes"
)
missing_manual <- setdiff(required_manual, names(manual_df))
if (length(missing_manual) > 0) {
  stop("Manual annotation file missing columns: ", paste(missing_manual, collapse = ", "))
}

# =========================
# Merge
# =========================
manifest <- manual_df %>%
  left_join(design_summary, by = "dataset_id") %>%
  left_join(rawmeta_summary, by = "dataset_id") %>%
  mutate(
    organism = dplyr::coalesce(na_if(organism, ""), organism_raw),
    tissue = dplyr::coalesce(tissue_design, tissue_raw),
    assay_type = dplyr::coalesce(na_if(assay_type, ""), assay_type_raw),
    time_hr = suppressWarnings(as.numeric(time_hr_design)),
    grouping_variable = dplyr::coalesce(grouping_variable, "group")
  ) %>%
  select(
    dataset_id,
    dataset_class,
    organism,
    tissue,
    assay_type,
    time_hr,
    treatment_type,
    grouping_variable,
    include_in_analysis,
    notes
  ) %>%
  arrange(dataset_id)

write_tsv(manifest, out_file)
log_message("Saved manifest: ", out_file)
log_message("Rows written: ", nrow(manifest))