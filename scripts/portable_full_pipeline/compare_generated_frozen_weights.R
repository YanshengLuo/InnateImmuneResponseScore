#!/usr/bin/env Rscript

this_file <- normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_file), "lib", "portable_pipeline_utils.R"))

ctx <- imrs_load_context(imrs_parse_args())
imrs_create_output_dirs(ctx)
log_msg <- imrs_log_function(ctx, "03_compare_generated_frozen_weights")
output_file <- file.path(ctx$paths$output_root, "manifests",
                         "frozen_gene_weight_comparison.tsv")
detail_file <- file.path(ctx$paths$output_root, "manifests",
                         "frozen_gene_weight_differences.tsv")
generated_file <- file.path(ctx$paths$output_root, "05_score", "anchors",
                            "gene_weights.tsv")
canonical_file <- file.path(ctx$paths$released_derived_root, "frozen_gene_weights.tsv")
tolerance <- as.numeric(ctx$config$comparison$released_weight_tolerance %||%
                          ctx$config$modeling$released_weight_tolerance %||% 1e-8)

if (!file.exists(generated_file)) {
  stop("Generated weights are missing: ", generated_file, call. = FALSE)
}

make_summary <- function(status, generated_count, canonical_count = NA_integer_,
                         overlapping_count = NA_integer_, missing_generated = NA_integer_,
                         missing_canonical = NA_integer_, correlation = NA_real_,
                         max_difference = NA_real_, mean_difference = NA_real_,
                         within = NA) {
  tibble::tibble(
    regenerated_gene_count = generated_count,
    canonical_gene_count = canonical_count,
    overlapping_gene_count = overlapping_count,
    missing_from_regenerated = missing_generated,
    missing_from_canonical = missing_canonical,
    coefficient_correlation = correlation,
    max_absolute_difference = max_difference,
    mean_absolute_difference = mean_difference,
    within_tolerance = within,
    status_message = status
  )
}

generated <- readr::read_tsv(generated_file, show_col_types = FALSE, progress = FALSE)
if (ctx$mode == "test") {
  summary <- make_summary(
    "TEST_ONLY_NON_CANONICAL: canonical equivalence comparison intentionally not performed.",
    nrow(generated)
  )
  readr::write_tsv(summary, output_file)
  log_msg(summary$status_message[[1]])
  quit(status = 0)
}
if (!file.exists(canonical_file)) {
  summary <- make_summary(
    "Canonical released frozen_gene_weights.tsv was not found; no equivalence claim made.",
    nrow(generated)
  )
  readr::write_tsv(summary, output_file)
  log_msg(summary$status_message[[1]])
  quit(status = 0)
}

canonical <- readr::read_tsv(canonical_file, show_col_types = FALSE, progress = FALSE)
prepare_weights <- function(df, label) {
  gene_col <- imrs_first_named_column(df, c("gene", "gene_id"))
  value_col <- imrs_first_named_column(df, c("weight", "beta_meta"))
  if (is.na(gene_col) || is.na(value_col)) {
    stop(label, " gene weights table must contain gene/gene_id and weight/beta_meta columns.",
         call. = FALSE)
  }
  tibble::tibble(
    gene = as.character(df[[gene_col]]),
    value = suppressWarnings(as.numeric(df[[value_col]]))
  )
}

regen <- prepare_weights(generated, "Generated") |>
  dplyr::rename(regenerated_coefficient = value)
canon <- prepare_weights(canonical, "Canonical") |>
  dplyr::rename(canonical_coefficient = value)
joined <- dplyr::full_join(canon, regen, by = "gene") |>
  dplyr::mutate(
    absolute_difference = abs(regenerated_coefficient - canonical_coefficient),
    comparison_status = dplyr::case_when(
      is.na(regenerated_coefficient) ~ "missing_from_regenerated",
      is.na(canonical_coefficient) ~ "missing_from_canonical",
      absolute_difference <= tolerance ~ "within_tolerance",
      TRUE ~ "different"
    )
  )
overlap <- joined |>
  dplyr::filter(!is.na(regenerated_coefficient), !is.na(canonical_coefficient))
correlation <- if (nrow(overlap) > 1) {
  stats::cor(overlap$regenerated_coefficient, overlap$canonical_coefficient,
             use = "complete.obs")
} else {
  NA_real_
}
within <- nrow(overlap) == nrow(joined) &&
  all(overlap$absolute_difference <= tolerance, na.rm = TRUE)
summary <- make_summary(
  if (within) "PASS: regenerated frozen weights match the released canonical table within tolerance."
  else "MISMATCH: regenerated frozen weights differ from the released canonical table; inspect the detail report.",
  nrow(regen),
  nrow(canon),
  nrow(overlap),
  sum(joined$comparison_status == "missing_from_regenerated"),
  sum(joined$comparison_status == "missing_from_canonical"),
  correlation,
  if (nrow(overlap) > 0) max(overlap$absolute_difference, na.rm = TRUE) else NA_real_,
  if (nrow(overlap) > 0) mean(overlap$absolute_difference, na.rm = TRUE) else NA_real_,
  within
)
readr::write_tsv(joined, detail_file)
readr::write_tsv(summary, output_file)
log_msg(summary$status_message[[1]])
