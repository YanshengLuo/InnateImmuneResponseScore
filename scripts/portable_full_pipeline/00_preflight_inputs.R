#!/usr/bin/env Rscript

this_file <- normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_file), "lib", "portable_pipeline_utils.R"))

args <- imrs_parse_args()
ctx <- imrs_load_context(args)
imrs_create_output_dirs(ctx)
log_msg <- imrs_log_function(ctx, "00_preflight")
out_dir <- file.path(ctx$paths$output_root, "00_preflight")
downstream_packages <- c("purrr", "tidyr", "DESeq2", "ggplot2", "pROC")
missing_downstream_packages <- downstream_packages[
  !vapply(downstream_packages, requireNamespace, logical(1), quietly = TRUE)
]

configured <- imrs_configured_datasets(ctx)
available <- imrs_available_datasets(ctx)
required <- if (ctx$mode == "canonical") {
  ctx$locked_anchors
} else if (ctx$mode == "all_scored") {
  ctx$locked_anchors
} else {
  available
}
label <- imrs_model_label(ctx)

log_msg("Preflight started in ", ctx$mode, " mode (", label, ").")
log_msg("Counts root: ", ctx$paths$counts_root)
log_msg("Verified metadata directory: ", ctx$paths$verified_metadata_dir)
log_msg("Configured datasets: ", paste(configured, collapse = ", "))

input_inventory <- dplyr::bind_rows(lapply(configured, function(dataset_id) {
  count_path <- imrs_count_path(ctx, dataset_id)
  metadata_path <- imrs_metadata_path(ctx, dataset_id)
  tibble::tibble(
    dataset_id = dataset_id,
    required_for_run = dataset_id %in% required,
    locked_anchor = dataset_id %in% ctx$locked_anchors,
    count_file = count_path,
    count_file_exists = file.exists(count_path),
    metadata_file = metadata_path,
    metadata_file_exists = file.exists(metadata_path),
    model_label = label
  )
}))

metadata_rows <- list()
sample_rows <- list()
issues <- character()
if (length(missing_downstream_packages) > 0) {
  issues <- c(issues, paste0(
    "Missing required downstream R packages: ",
    paste(missing_downstream_packages, collapse = ", ")
  ))
}
for (dataset_id in configured) {
  row <- input_inventory[input_inventory$dataset_id == dataset_id, , drop = FALSE]
  if (!row$metadata_file_exists) {
    metadata_rows[[length(metadata_rows) + 1L]] <- tibble::tibble(
      dataset_id = dataset_id, metadata_file = row$metadata_file,
      metadata_file_exists = FALSE, n_samples = NA_integer_,
      n_control = NA_integer_, n_delivery = NA_integer_,
      labels_generated = FALSE, issue = "metadata_file_missing",
      model_label = label
    )
    if (dataset_id %in% required) issues <- c(issues, paste0(dataset_id, ": metadata file missing"))
    next
  }
  design <- tryCatch(imrs_scoring_design(ctx, dataset_id), error = function(e) e)
  if (inherits(design, "error")) {
    issue <- conditionMessage(design)
    metadata_rows[[length(metadata_rows) + 1L]] <- tibble::tibble(
      dataset_id = dataset_id, metadata_file = row$metadata_file,
      metadata_file_exists = TRUE, n_samples = NA_integer_,
      n_control = NA_integer_, n_delivery = NA_integer_,
      labels_generated = FALSE, issue = issue, model_label = label
    )
    if (dataset_id %in% required) issues <- c(issues, paste0(dataset_id, ": ", issue))
    next
  }
  metadata_rows[[length(metadata_rows) + 1L]] <- tibble::tibble(
    dataset_id = dataset_id, metadata_file = row$metadata_file,
    metadata_file_exists = TRUE, n_samples = nrow(design),
    n_control = sum(design$condition_simple == "CONTROL"),
    n_delivery = sum(design$condition_simple == "DELIVERY"),
    labels_generated = all(design$condition_simple %in% c("CONTROL", "DELIVERY")),
    issue = NA_character_, model_label = label
  )
  if (!row$count_file_exists) {
    sample_rows[[length(sample_rows) + 1L]] <- tibble::tibble(
      dataset_id = dataset_id, count_file = row$count_file, n_metadata_samples = nrow(design),
      n_count_sample_columns = NA_integer_, n_overlapping_samples = 0L,
      gene_identifier_column = NA_character_, raw_integer_counts = FALSE,
      match_status = "count_file_missing", model_label = label
    )
    if (dataset_id %in% required) issues <- c(issues, paste0(dataset_id, ": count file missing"))
    next
  }
  count_check <- tryCatch({
    obj <- imrs_read_counts(row$count_file, design$sample_id, strip_versions = FALSE)
    imrs_assert_raw_counts(obj$matrix, row$count_file)
    obj
  }, error = function(e) e)
  if (inherits(count_check, "error")) {
    issue <- conditionMessage(count_check)
    sample_rows[[length(sample_rows) + 1L]] <- tibble::tibble(
      dataset_id = dataset_id, count_file = row$count_file, n_metadata_samples = nrow(design),
      n_count_sample_columns = NA_integer_, n_overlapping_samples = NA_integer_,
      gene_identifier_column = NA_character_, raw_integer_counts = FALSE,
      match_status = issue, model_label = label
    )
    if (dataset_id %in% required) issues <- c(issues, paste0(dataset_id, ": ", issue))
  } else {
    overlap <- length(count_check$sample_columns)
    status <- if (overlap == nrow(design)) "all_metadata_samples_matched" else
      "partial_metadata_sample_overlap"
    sample_rows[[length(sample_rows) + 1L]] <- tibble::tibble(
      dataset_id = dataset_id, count_file = row$count_file, n_metadata_samples = nrow(design),
      n_count_sample_columns = ncol(count_check$matrix), n_overlapping_samples = overlap,
      gene_identifier_column = count_check$gene_column, raw_integer_counts = TRUE,
      match_status = status, model_label = label
    )
    if (dataset_id %in% required && overlap == 0) {
      issues <- c(issues, paste0(dataset_id, ": no overlapping samples"))
    }
  }
}

metadata_inventory <- dplyr::bind_rows(metadata_rows)
sample_report <- dplyr::bind_rows(sample_rows)
readr::write_tsv(input_inventory, file.path(out_dir, "input_inventory.tsv"))
readr::write_tsv(metadata_inventory, file.path(out_dir, "metadata_inventory.tsv"))
readr::write_tsv(sample_report, file.path(out_dir, "sample_match_report.tsv"))

if (ctx$mode %in% c("canonical", "all_scored")) {
  missing_anchors <- setdiff(ctx$locked_anchors, available)
  if (length(missing_anchors) > 0) {
    issues <- unique(c(issues, paste0(
      ctx$mode, " mode lacks required locked anchors: ",
      paste(missing_anchors, collapse = ", ")
    )))
  }
  if (ctx$mode == "all_scored") {
    log_msg("ALL_SCORED MODE: frozen weights are reconstructed from locked anchors only; ",
            "all other available datasets are scored without refitting.")
  }
} else {
  log_msg("TEST MODE: outputs are non-canonical/test-only and cannot establish canonical equivalence.")
  if (length(available) == 0) issues <- c(issues, "Test mode has no complete dataset inputs.")
}

passed <- length(unique(issues)) == 0
status_lines <- c(
  paste0("mode: ", ctx$mode),
  paste0("model_label: ", label),
  paste0("status: ", if (passed) "PASS" else "FAIL"),
  paste0("available_datasets: ", paste(available, collapse = ", ")),
  paste0("locked_anchors_required: ", paste(ctx$locked_anchors, collapse = ", ")),
  paste0("required_downstream_packages: ", paste(downstream_packages, collapse = ", ")),
  paste0("missing_downstream_packages: ",
         if (length(missing_downstream_packages) == 0) "none" else
           paste(missing_downstream_packages, collapse = ", ")),
  if (ctx$mode == "test") "canonical_equivalence_claim: prohibited_test_only" else
    "canonical_equivalence_claim: pending_locked_anchor_reconstructed_weight_comparison",
  if (length(issues) > 0) paste0("issues: ", paste(unique(issues), collapse = " | ")) else
    "issues: none"
)
writeLines(status_lines, file.path(out_dir, "preflight_status.txt"), useBytes = TRUE)
log_msg("Preflight status: ", if (passed) "PASS" else "FAIL")

if (!passed && !ctx$dry_run) {
  stop("Preflight failed: ", paste(unique(issues), collapse = " | "), call. = FALSE)
}
if (!passed && ctx$dry_run) {
  log_msg("Dry-run preflight recorded missing requirements without starting downstream stages.")
} else {
  log_msg("Preflight completed successfully for datasets: ",
          paste(imrs_active_datasets(ctx), collapse = ", "))
}
