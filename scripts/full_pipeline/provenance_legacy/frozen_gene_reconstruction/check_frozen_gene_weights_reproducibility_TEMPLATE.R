#!/usr/bin/env Rscript

# Optional template: compare regenerated frozen gene weights against the
# released frozen_gene_weights.tsv. This script is intentionally not called by
# run_all_manuscript_outputs_v6.R.

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

read_simple_yaml <- function(file) {
  if (!file.exists(file)) return(list())
  lines <- readLines(file, warn = FALSE)
  lines <- lines[grepl(":", lines) & !grepl("^\\s*#", lines)]
  out <- list()
  for (line in lines) {
    key <- trimws(sub(":.*$", "", line))
    value <- trimws(sub("^[^:]+:", "", line))
    value <- trimws(gsub("^['\"]|['\"]$", "", value))
    out[[key]] <- value
  }
  out
}

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                      error = function(e) NA_character_)
repo_root <- if (!is.na(this_file)) {
  normalizePath(file.path(dirname(this_file), "..", "..", ".."),
                winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
if (!dir.exists(file.path(repo_root, "config"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

config_file <- file.path(repo_root, "config", "config.yml")
if (!file.exists(config_file)) {
  config_file <- file.path(repo_root, "config", "config_template.yml")
}
config <- read_simple_yaml(config_file)

resolve_path <- function(value, default) {
  if (is.null(value) || !nzchar(value)) value <- default
  if (grepl("^([A-Za-z]:|/)", value)) {
    normalizePath(value, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(repo_root, value), winslash = "/", mustWork = FALSE)
  }
}

released_file <- resolve_path(config$frozen_gene_weights %||% "data_release_templates/derived/frozen_gene_weights.tsv",
                              "data_release_templates/derived/frozen_gene_weights.tsv")
regenerated_file <- resolve_path(config$regenerated_frozen_gene_weights %||% "results_release_templates/manifests/regenerated_gene_weights.tsv",
                                 "results_release_templates/manifests/regenerated_gene_weights.tsv")
out_file <- resolve_path("results_release_templates/manifests/frozen_gene_weights_reproducibility_check.tsv",
                         "results_release_templates/manifests/frozen_gene_weights_reproducibility_check.tsv")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write_result <- function(row) {
  utils::write.table(row, out_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

if (!file.exists(released_file)) {
  write_result(data.frame(
    released_file = released_file,
    regenerated_file = regenerated_file,
    n_genes_released = NA_integer_,
    n_genes_regenerated = NA_integer_,
    n_gene_ids_matching = NA_integer_,
    all_gene_ids_match_yes_no = "no",
    all_weights_identical_yes_no = "no",
    max_abs_weight_difference = NA_real_,
    coefficient_column_compared = NA_character_,
    status = "released_file_missing",
    notes = "Stage data_release_templates/derived/frozen_gene_weights.tsv before running this optional check."
  ))
  quit(status = 0)
}

if (!file.exists(regenerated_file)) {
  released <- utils::read.delim(released_file, check.names = FALSE)
  write_result(data.frame(
    released_file = released_file,
    regenerated_file = regenerated_file,
    n_genes_released = nrow(released),
    n_genes_regenerated = NA_integer_,
    n_gene_ids_matching = NA_integer_,
    all_gene_ids_match_yes_no = "not_run",
    all_weights_identical_yes_no = "not_run",
    max_abs_weight_difference = NA_real_,
    coefficient_column_compared = NA_character_,
    status = "not_run",
    notes = "No regenerated frozen weights file was found. Run the optional full reconstruction path, then set regenerated_frozen_gene_weights in config.yml or place the file at the default path."
  ))
  quit(status = 0)
}

released <- utils::read.delim(released_file, check.names = FALSE)
regenerated <- utils::read.delim(regenerated_file, check.names = FALSE)

gene_col <- intersect(c("gene", "gene_id", "ensembl_gene_id"), names(released))[1]
regen_gene_col <- intersect(c(gene_col, "gene", "gene_id", "ensembl_gene_id"), names(regenerated))[1]
coef_col <- intersect(c("weight", "beta_meta", "coefficient"), names(released))[1]
regen_coef_col <- intersect(c(coef_col, "weight", "beta_meta", "coefficient"), names(regenerated))[1]

if (is.na(gene_col) || is.na(regen_gene_col) || is.na(coef_col) || is.na(regen_coef_col)) {
  write_result(data.frame(
    released_file = released_file,
    regenerated_file = regenerated_file,
    n_genes_released = nrow(released),
    n_genes_regenerated = nrow(regenerated),
    n_gene_ids_matching = NA_integer_,
    all_gene_ids_match_yes_no = "no",
    all_weights_identical_yes_no = "no",
    max_abs_weight_difference = NA_real_,
    coefficient_column_compared = paste(coef_col, regen_coef_col, sep = " vs "),
    status = "column_mismatch",
    notes = "Could not identify comparable gene or coefficient columns."
  ))
  quit(status = 0)
}

rel <- released[, c(gene_col, coef_col)]
reg <- regenerated[, c(regen_gene_col, regen_coef_col)]
names(rel) <- c("gene_id", "released_weight")
names(reg) <- c("gene_id", "regenerated_weight")
merged <- merge(rel, reg, by = "gene_id", all = TRUE)
diffs <- abs(as.numeric(merged$released_weight) - as.numeric(merged$regenerated_weight))

write_result(data.frame(
  released_file = released_file,
  regenerated_file = regenerated_file,
  n_genes_released = nrow(released),
  n_genes_regenerated = nrow(regenerated),
  n_gene_ids_matching = sum(!is.na(merged$released_weight) & !is.na(merged$regenerated_weight)),
  all_gene_ids_match_yes_no = ifelse(setequal(rel$gene_id, reg$gene_id), "yes", "no"),
  all_weights_identical_yes_no = ifelse(all(diffs == 0, na.rm = TRUE) && !any(is.na(diffs)), "yes", "no"),
  max_abs_weight_difference = ifelse(length(diffs) == 0 || all(is.na(diffs)), NA_real_, max(diffs, na.rm = TRUE)),
  coefficient_column_compared = paste(coef_col, regen_coef_col, sep = " vs "),
  status = "completed",
  notes = "Optional comparison only; this script does not refit IMRS and is not part of the default reviewer manuscript runner."
))
