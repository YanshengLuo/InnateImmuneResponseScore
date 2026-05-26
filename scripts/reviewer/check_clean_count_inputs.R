#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0L) {
  stop("This script must be run with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]),
                             winslash = "/", mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                           winslash = "/", mustWork = FALSE)

datasets <- c(
  "GSE119119",
  "GSE139529",
  "GSE166655",
  "GSE167521",
  "GSE178313",
  "GSE262515",
  "GSE264344",
  "GSE279372",
  "GSE279743",
  "GSE279744",
  "GSE314070",
  "GSE39129"
)

inspect_count_file <- function(dataset_id) {
  path <- file.path(repo_root, "data", "counts", dataset_id, "featurecounts",
                    "validation", "gene_counts_clean.tsv")
  exists <- file.exists(path)
  bytes <- if (exists) file.info(path)$size else NA_real_
  non_empty <- isTRUE(exists) && !is.na(bytes) && bytes > 0
  header_ok <- FALSE
  status <- "MISSING"

  if (non_empty) {
    header <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                       error = function(e) character())
    header_ok <- length(header) == 1L && nzchar(header[[1]]) &&
      grepl("\t", header[[1]], fixed = TRUE)
    status <- if (header_ok) "PASS" else "FAIL_HEADER"
  } else if (exists) {
    status <- "FAIL_EMPTY"
  }

  data.frame(
    dataset_id = dataset_id,
    path = path,
    file_size_mb = if (is.na(bytes)) NA_real_ else round(bytes / (1024^2), 2),
    exists = exists,
    non_empty = non_empty,
    header_readable_tsv = header_ok,
    status = status,
    stringsAsFactors = FALSE
  )
}

report <- do.call(rbind, lapply(datasets, inspect_count_file))
preferred_output_dir <- file.path(repo_root, "results", "full_pipeline_v6", "00_preflight")
output_dir <- if (dir.exists(preferred_output_dir)) {
  preferred_output_dir
} else {
  file.path(repo_root, "results")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "clean_count_input_check.tsv")

write.table(report, output_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
print(report, row.names = FALSE)
cat("\nWrote clean-count input report to: ", output_path, "\n", sep = "")

if (!all(report$status == "PASS")) {
  stop("One or more expected clean count inputs failed validation.", call. = FALSE)
}
