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

counts_root <- file.path(imrs_root, "03_counts")
sample_metadata_file <- file.path(framework_root, "inputs", "02_sample_metadata.tsv")

out_index_file <- file.path(framework_root, "inputs", "03_expression_matrix_index.tsv")
out_matrix_dir <- file.path(framework_root, "inputs", "expression_matrices")
log_file <- file.path(framework_root, "logs", "03_build_expression_matrix_index.log")

dir.create(dirname(out_index_file), recursive = TRUE, showWarnings = FALSE)
dir.create(out_matrix_dir, recursive = TRUE, showWarnings = FALSE)
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
read_counts_matrix <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (ncol(df) < 2) stop("Counts file has fewer than 2 columns: ", path)
  names(df)[1] <- "gene_id"
  df
}

# =========================
# Inputs
# =========================
if (!file.exists(sample_metadata_file)) {
  stop("Missing sample metadata file: ", sample_metadata_file)
}

sample_metadata <- read_tsv(sample_metadata_file, show_col_types = FALSE, progress = FALSE)

required_cols <- c("dataset_id", "sample_id")
missing_cols <- setdiff(required_cols, names(sample_metadata))
if (length(missing_cols) > 0) {
  stop("02_sample_metadata.tsv missing columns: ", paste(missing_cols, collapse = ", "))
}

dataset_ids <- sample_metadata %>%
  distinct(dataset_id) %>%
  pull(dataset_id)

log_message("Counts root: ", counts_root)
log_message("Datasets in sample metadata: ", length(dataset_ids))

index_rows <- list()

# =========================
# Main loop
# =========================
for (dataset_id in dataset_ids) {
  counts_file <- file.path(
    counts_root,
    dataset_id,
    "featurecounts",
    "validation",
    "gene_counts_clean.tsv"
  )

  if (!file.exists(counts_file)) {
    log_message("Missing counts file for ", dataset_id, ": ", counts_file)
    next
  }

  log_message("Reading counts for ", dataset_id, ": ", counts_file)

  counts_df <- tryCatch(
    read_counts_matrix(counts_file),
    error = function(e) e
  )

  if (inherits(counts_df, "error")) {
    log_message("Failed to read counts for ", dataset_id, ": ", conditionMessage(counts_df))
    next
  }

  expected_samples <- sample_metadata %>%
    filter(dataset_id == !!dataset_id) %>%
    pull(sample_id) %>%
    unique()

  available_samples <- intersect(expected_samples, colnames(counts_df))

  if (length(available_samples) == 0) {
    log_message("No matching sample columns found for ", dataset_id)
    next
  }

  missing_samples <- setdiff(expected_samples, colnames(counts_df))
  if (length(missing_samples) > 0) {
    log_message("Missing sample columns in counts for ", dataset_id, ": ",
                paste(missing_samples, collapse = ", "))
  }

  ordered_samples <- expected_samples[expected_samples %in% available_samples]

  out_df <- counts_df %>%
    select(gene_id, all_of(ordered_samples))

  out_file <- file.path(out_matrix_dir, paste0(dataset_id, "_expression_matrix.tsv"))
  write_tsv(out_df, out_file)

  index_rows[[length(index_rows) + 1]] <- tibble(
    dataset_id = dataset_id,
    expression_matrix_file = file.path("inputs", "expression_matrices", paste0(dataset_id, "_expression_matrix.tsv")),
    value_type = "raw_counts",
    gene_id_type = "as_is",
    n_genes = nrow(out_df),
    n_samples = ncol(out_df) - 1
  )

  log_message("Saved standardized matrix for ", dataset_id, " with ",
              nrow(out_df), " genes and ", ncol(out_df) - 1, " samples")
}

index_df <- bind_rows(index_rows)

if (nrow(index_df) == 0) {
  stop("No expression matrices were successfully built.")
}

write_tsv(index_df, out_index_file)
log_message("Saved expression matrix index: ", out_index_file)
log_message("Datasets written: ", nrow(index_df))