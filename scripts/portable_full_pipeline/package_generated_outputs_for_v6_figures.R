#!/usr/bin/env Rscript

this_file <- normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_file), "lib", "portable_pipeline_utils.R"))

imrs_require(c("yaml", "readr", "dplyr", "tibble"))
ctx <- imrs_load_context(imrs_parse_args())
imrs_create_output_dirs(ctx)
log_msg <- imrs_log_function(ctx, "06_package_v6_figure_inputs")
label <- imrs_model_label(ctx)
generated_root <- ctx$paths$generated_derived_root
released_root <- ctx$paths$released_derived_root
if (identical(imrs_normalize(generated_root, NULL, FALSE),
              imrs_normalize(released_root, NULL, FALSE))) {
  stop("Generated derived root resolves to released data/derived. Refusing to overwrite released inputs.",
       call. = FALSE)
}
figure_root <- file.path(generated_root, "figure_inputs")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
manifest_rows <- list()

stage_file <- function(source, destination, action, required = TRUE) {
  if (!file.exists(source)) {
    if (required) stop("Required figure-packaging input missing: ", source, call. = FALSE)
    manifest_rows[[length(manifest_rows) + 1L]] <<- tibble::tibble(
      source_file = source, packaged_file = destination, action = action,
      status = "source_not_available", model_label = label
    )
    return(invisible(FALSE))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(source, destination, overwrite = TRUE)
  if (!ok) stop("Could not write staged derived input: ", destination, call. = FALSE)
  manifest_rows[[length(manifest_rows) + 1L]] <<- tibble::tibble(
    source_file = source, packaged_file = destination, action = action,
    status = "written", model_label = label
  )
  invisible(TRUE)
}

anchors <- file.path(ctx$paths$output_root, "05_score", "anchors")
evaluation <- file.path(ctx$paths$output_root, "05_score", "transfer", "eval")
regenerated <- stats::setNames(
  c(
    file.path(generated_root, "frozen_gene_weights.tsv"),
    file.path(generated_root, "gene_power.tsv"),
    file.path(figure_root, "gene_weights.tsv"),
    file.path(figure_root, "gene_heterogeneity.tsv"),
    file.path(figure_root, "gene_power.tsv"),
    file.path(figure_root, "support_by_dataset.tsv"),
    file.path(figure_root, "contrast_counts_by_dataset.tsv"),
    file.path(figure_root, "contrast_counts_by_dataset_6B.tsv"),
    file.path(figure_root, "step09_split_eval.tsv"),
    file.path(figure_root, "step09_split_summary.tsv"),
    file.path(figure_root, "step09_split_sample_level.tsv")
  ),
  c(
    file.path(anchors, "gene_weights.tsv"),
    file.path(anchors, "gene_power.tsv"),
    file.path(anchors, "gene_weights.tsv"),
    file.path(anchors, "gene_heterogeneity.tsv"),
    file.path(anchors, "gene_power.tsv"),
    file.path(anchors, "support_by_dataset.tsv"),
    file.path(anchors, "contrast_counts_by_dataset.tsv"),
    file.path(anchors, "contrast_counts_by_dataset_6B.tsv"),
    file.path(evaluation, "step09_split_eval.tsv"),
    file.path(evaluation, "step09_split_summary.tsv"),
    file.path(evaluation, "step09_split_sample_level.tsv")
  )
)
for (i in seq_along(regenerated)) {
  stage_file(names(regenerated)[[i]], regenerated[[i]], "regenerated_count_level_output")
}

mapping_source <- file.path(released_root, "gene_symbol_mapping.tsv")
stage_file(mapping_source, file.path(generated_root, "gene_symbol_mapping.tsv"),
           "released_mapping_copy_through", required = FALSE)
stage_file(mapping_source, file.path(figure_root, "gene_symbol_mapping.tsv"),
           "released_mapping_copy_through", required = FALSE)

release_figure_root <- file.path(released_root, "figure_inputs")
generated_names <- basename(unname(regenerated[dirname(unname(regenerated)) == figure_root]))
if (dir.exists(release_figure_root)) {
  released_figure_files <- list.files(release_figure_root, pattern = "\\.tsv$",
                                      full.names = TRUE)
  released_figure_files <- released_figure_files[
    !basename(released_figure_files) %in% generated_names &
      basename(released_figure_files) != "gene_symbol_mapping.tsv"
  ]
  for (source in released_figure_files) {
    stage_file(source, file.path(figure_root, basename(source)),
               "released_robustness_or_context_copy_through", required = FALSE)
  }
}

released_top_level <- c(
  "supplement_dataset_split_provenance_v7.tsv",
  "weak_dataset_paper_context_audit.tsv",
  "label_permutation_null_summary.tsv",
  "leave_one_gene_out_summary.tsv",
  "gene_dominance_summary.tsv",
  "threshold_sensitivity_summary.tsv",
  "leave_one_anchor_out_summary.tsv",
  "baseline_signature_summary_by_group.tsv",
  "coefficient_sensitivity_summary.tsv"
)
for (filename in released_top_level) {
  stage_file(file.path(released_root, filename), file.path(generated_root, filename),
             "released_robustness_or_context_copy_through", required = FALSE)
}

manifest <- dplyr::bind_rows(manifest_rows) |>
  dplyr::mutate(
    canonical_released_tables_overwritten = FALSE,
    packaged_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
manifest_path <- file.path(ctx$paths$output_root, "manifests",
                           "v6_figure_input_packaging_manifest.tsv")
readr::write_tsv(manifest, manifest_path)

manuscript_output_root <- if (ctx$mode == "test") {
  file.path(ctx$paths$output_root, "manuscript_outputs")
} else {
  imrs_normalize(
    ctx$config$manuscript$output_root %||% file.path(ctx$paths$output_root, "manuscript_outputs"),
    ctx$project_root, FALSE
  )
}
staged_config <- list(
  project_root = ctx$repo_root,
  figure_input_dir = figure_root,
  derived_input_dir = generated_root,
  derived_data_dir = generated_root,
  results_dir = manuscript_output_root,
  logs_dir = file.path(manuscript_output_root, "logs"),
  manifests_dir = file.path(manuscript_output_root, "manifests"),
  figures_dir = file.path(manuscript_output_root, "figures"),
  supplementary_figures_dir = file.path(manuscript_output_root, "supplementary_figures"),
  supplementary_tables_dir = file.path(manuscript_output_root, "supplementary_tables"),
  manuscript_output_dir = manuscript_output_root,
  priority3_enrichment_dir = file.path(manuscript_output_root, "priority3_gene_program_enrichment"),
  execute_active_scripts = TRUE,
  run_internal_readiness = FALSE,
  figure_helper_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                   "figure_helpers_v6.R"),
  panel_builder_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                   "panel_builders_v6.R"),
  workflow_diagram_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                      "merged_workflow_v6.R"),
  frozen_gene_weights = file.path(generated_root, "frozen_gene_weights.tsv"),
  gene_symbol_mapping = file.path(generated_root, "gene_symbol_mapping.tsv"),
  gene_power = file.path(generated_root, "gene_power.tsv"),
  gene_heterogeneity = file.path(figure_root, "gene_heterogeneity.tsv"),
  support_by_dataset = file.path(figure_root, "support_by_dataset.tsv"),
  core_support_all_genes = file.path(figure_root, "support_by_dataset.tsv"),
  provenance_table = file.path(generated_root, "supplement_dataset_split_provenance_v7.tsv"),
  manuscript_role_table = file.path(figure_root, "manuscript_dataset_role_table.tsv"),
  boundary_audit_table = file.path(generated_root, "weak_dataset_paper_context_audit.tsv"),
  label_permutation_summary = file.path(generated_root, "label_permutation_null_summary.tsv"),
  leave_one_gene_out_summary = file.path(generated_root, "leave_one_gene_out_summary.tsv"),
  gene_dominance_summary = file.path(generated_root, "gene_dominance_summary.tsv"),
  threshold_sensitivity_summary = file.path(generated_root, "threshold_sensitivity_summary.tsv"),
  leave_one_anchor_out_summary = file.path(generated_root, "leave_one_anchor_out_summary.tsv"),
  comparator_benchmarking_summary = file.path(generated_root, "baseline_signature_summary_by_group.tsv"),
  coefficient_sensitivity_summary = file.path(generated_root, "coefficient_sensitivity_summary.tsv"),
  step09_split_eval = file.path(figure_root, "step09_split_eval.tsv"),
  step09_split_summary = file.path(figure_root, "step09_split_summary.tsv"),
  source_v5_dir = generated_root
)
staged_config_path <- file.path(ctx$paths$output_root, "manifests",
                                "manuscript_generated_inputs_config.yml")
yaml::write_yaml(staged_config, staged_config_path)
log_msg("Packaged regenerated figure inputs under: ", generated_root)
log_msg("Released data/derived tables were not overwritten.")
log_msg("Staged v6 manuscript config: ", staged_config_path)
