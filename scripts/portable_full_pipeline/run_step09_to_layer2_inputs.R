#!/usr/bin/env Rscript

this_file <- normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)
script_root <- dirname(this_file)
source(file.path(script_root, "lib", "portable_pipeline_utils.R"))

args <- imrs_parse_args()
if (isTRUE(args$help)) {
  cat(
    "Usage: Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R",
    " --config config/full_pipeline_config.yml --mode all_scored [--dry-run] [--force]\n"
  )
  quit(status = 0)
}
ctx <- imrs_load_context(args)
if (ctx$mode != "all_scored") {
  stop("The Step09-to-Layer2 bridge requires --mode all_scored because it needs ",
       "validation score/evaluation tables generated with frozen locked-anchor weights.",
       call. = FALSE)
}
imrs_require(c("yaml", "readr", "dplyr", "tibble"))

output_root <- ctx$paths$output_root
bridge_root <- file.path(output_root, "layer2_generated_inputs")
derived_root <- file.path(bridge_root, "data", "derived")
figure_root <- file.path(derived_root, "figure_inputs")
work_root <- file.path(bridge_root, "work")
audit_root <- file.path(work_root, "audit")
audit_results <- file.path(audit_root, "results")
audit_report <- file.path(audit_root, "report")
extra_root <- file.path(work_root, "publication_extra_analyses")
metadata_stage <- file.path(work_root, "verified_metadata")
logs_dir <- file.path(bridge_root, "logs")
manifest_path <- file.path(bridge_root, "layer2_bridge_manifest.tsv")
layer2_config_path <- file.path(bridge_root, "layer2_generated_inputs_config.yml")
ported_bridge <- file.path(script_root, "bridge_to_layer2")

required_generated <- c(
  file.path(output_root, "05_score", "anchors", "gene_weights.tsv"),
  file.path(output_root, "05_score", "anchors", "gene_power.tsv"),
  file.path(output_root, "05_score", "anchors", "gene_heterogeneity.tsv"),
  file.path(output_root, "05_score", "anchors", "support_by_dataset.tsv"),
  file.path(output_root, "05_score", "transfer", "eval", "step09_split_eval.tsv"),
  file.path(output_root, "05_score", "transfer", "eval", "step09_split_summary.tsv"),
  file.path(output_root, "05_score", "transfer", "eval", "step09_split_sample_level.tsv"),
  file.path(output_root, "01_designs", "scoring"),
  file.path(output_root, "01_designs", "splited"),
  file.path(output_root, "04_de", "comparison")
)
missing_generated <- required_generated[!file.exists(required_generated) &
                                          !dir.exists(required_generated)]
if (length(missing_generated) > 0) {
  stop("Bridge prerequisites are missing. Run the count-level pipeline with --mode all_scored ",
       "through Step09 before running this bridge. Missing: ",
       paste(missing_generated, collapse = "; "), call. = FALSE)
}

audit_scripts <- file.path(ported_bridge, "audit", c(
  "01_dataset_audit_table.R",
  "03_leave_one_gene_out.R",
  "04_gene_dominance.R",
  "05_threshold_sensitivity.R",
  "06_leave_one_anchor_out.R",
  "07_weak_dataset_paper_context_audit.R",
  "08_manuscript_ready_cleanup.R",
  "14_create_supplement_provenance_v7.R"
))
extra_scripts <- file.path(ported_bridge, "publication_extra", c(
  "01_label_permutation_null.R",
  "02_baseline_signature_benchmarking.R",
  "03_coefficient_sensitivity_summary.R"
))
scripts <- c(audit_scripts, extra_scripts)
missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts) > 0) {
  stop("Ported bridge script missing: ", paste(missing_scripts, collapse = "; "),
       call. = FALSE)
}

cat("IMRS Step09-to-Layer2 bridge plan\n")
cat("  input Step09 root: ", file.path(output_root, "05_score", "transfer", "eval"), "\n",
    sep = "")
cat("  output package root: ", bridge_root, "\n", sep = "")
cat("  generated tables will not overwrite: ", ctx$paths$released_derived_root, "\n", sep = "")
cat("  copied original execution order:\n")
cat(paste0("    - ", vapply(scripts, imrs_rel_or_path, character(1), root = ctx$repo_root)),
    sep = "\n")
cat("\n")
if (ctx$dry_run) {
  cat("Dry run complete. No bridge analyses or output packaging executed.\n")
  quit(status = 0)
}
if (dir.exists(bridge_root) && file.exists(manifest_path) && !isTRUE(args$force)) {
  stop("A generated Layer2 bridge package already exists at ", bridge_root,
       ". Use --force to replace generated bridge products.", call. = FALSE)
}
if (identical(imrs_normalize(derived_root, NULL, FALSE),
              imrs_normalize(ctx$paths$released_derived_root, NULL, FALSE))) {
  stop("Generated Layer2 derived root resolves to committed released data/derived. Refusing to write.",
       call. = FALSE)
}

dirs <- c(figure_root, audit_results, audit_report, file.path(audit_root, "datasets"),
          extra_root, metadata_stage, logs_dir)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
log_file <- file.path(logs_dir, "run_step09_to_layer2_inputs.log")
writeLines(character(), log_file)
log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 paste(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

copy_tree <- function(source, destination) {
  if (!dir.exists(source)) stop("Directory required for bridge staging is missing: ", source,
                                call. = FALSE)
  files <- list.files(source, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[!dir.exists(files)]
  if (length(files) == 0) return(invisible(character()))
  source_norm <- paste0(imrs_normalize(source, NULL, FALSE), "/")
  normalized_files <- vapply(files, imrs_normalize, character(1), base = NULL,
                             must_work = FALSE)
  rel <- substring(normalized_files, nchar(source_norm) + 1L)
  targets <- file.path(destination, rel)
  invisible(lapply(dirname(targets), dir.create, recursive = TRUE, showWarnings = FALSE))
  ok <- file.copy(files, targets, overwrite = TRUE)
  if (!all(ok)) stop("Could not stage bridge input directory: ", source, call. = FALSE)
  targets
}

metadata_files <- list.files(ctx$paths$verified_metadata_dir, pattern = "\\.tsv$",
                             full.names = TRUE)
if (length(metadata_files) == 0) {
  stop("No verified metadata TSV files found for bridge provenance staging: ",
       ctx$paths$verified_metadata_dir, call. = FALSE)
}
if (!all(file.copy(metadata_files, file.path(metadata_stage, basename(metadata_files)),
                   overwrite = TRUE))) {
  stop("Could not stage verified metadata for supplement provenance generation.", call. = FALSE)
}
invisible(copy_tree(file.path(output_root, "01_designs", "scoring"),
                    file.path(metadata_stage, "scoring")))
invisible(copy_tree(file.path(output_root, "01_designs", "splited"),
                    file.path(metadata_stage, "splited")))
mapping_source <- file.path(ported_bridge, "reference_inputs", "weak_dataset_paper_mappings.tsv")
if (!file.copy(mapping_source, file.path(audit_root, "datasets",
                                         "weak_dataset_paper_mappings.tsv"), overwrite = TRUE)) {
  stop("Could not stage weak-dataset paper mapping reference for bridge execution.",
       call. = FALSE)
}

rscript <- if (.Platform$OS.type == "windows") file.path(R.home("bin"), "Rscript.exe") else
  file.path(R.home("bin"), "Rscript")
bridge_env <- c(
  IMRS_PROJECT_ROOT = output_root,
  IMRS_PORTED_OUTPUT_ROOT = output_root,
  IMRS_PORTED_COUNTS_ROOT = ctx$paths$counts_root,
  IMRS_PORTED_SCORING_DESIGN_ROOT = file.path(output_root, "01_designs", "scoring"),
  IMRS_PORTED_SPLIT_DESIGN_ROOT = file.path(output_root, "01_designs", "splited"),
  IMRS_PORTED_GENE_SYMBOL_MAP = file.path(ctx$paths$released_derived_root,
                                         "gene_symbol_mapping.tsv"),
  IMRS_LAYER2_AUDIT_ROOT = audit_root,
  IMRS_LAYER2_AUDIT_RESULTS_DIR = audit_results,
  IMRS_LAYER2_AUDIT_REPORT_DIR = audit_report,
  IMRS_LAYER2_EXTRA_ROOT = extra_root,
  IMRS_BRIDGE_VERIFIED_METADATA_ROOT = metadata_stage
)
prior <- Sys.getenv(names(bridge_env), unset = NA_character_)
on.exit({
  present <- !is.na(prior)
  if (any(present)) do.call(Sys.setenv, as.list(prior[present]))
  if (any(!present)) Sys.unsetenv(names(prior)[!present])
}, add = TRUE)
do.call(Sys.setenv, as.list(bridge_env))

run_one <- function(script) {
  log_msg("Running copied original bridge script: ", imrs_rel_or_path(script, ctx$repo_root))
  output <- system2(rscript, args = shQuote(script), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (length(output) > 0) cat(paste(output, collapse = "\n"), "\n", file = log_file, append = TRUE)
  if (status != 0L) {
    stop("Bridge script failed: ", script, ". Review ", log_file, ".", call. = FALSE)
  }
}
invisible(lapply(scripts, run_one))

manifest_rows <- list()
package_file <- function(source, destination, generated_or_reference, notes) {
  if (!file.exists(source)) stop("Required Layer2 packaging source missing: ", source, call. = FALSE)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE)) {
    stop("Could not package Layer2 input: ", destination, call. = FALSE)
  }
  tab <- tryCatch(readr::read_tsv(destination, show_col_types = FALSE, progress = FALSE),
                  error = function(e) NULL)
  manifest_rows[[length(manifest_rows) + 1L]] <<- tibble::tibble(
    output_file = imrs_rel_or_path(destination, bridge_root),
    source_files = source,
    generated_or_reference = generated_or_reference,
    row_count = if (is.null(tab)) NA_integer_ else nrow(tab),
    column_count = if (is.null(tab)) NA_integer_ else ncol(tab),
    checksum_or_size = unname(tools::md5sum(destination)),
    notes = notes
  )
}

anchors <- file.path(output_root, "05_score", "anchors")
evaluation <- file.path(output_root, "05_score", "transfer", "eval")
generated_top <- c(
  frozen_gene_weights.tsv = file.path(anchors, "gene_weights.tsv"),
  gene_power.tsv = file.path(anchors, "gene_power.tsv"),
  supplement_dataset_split_provenance_v7.tsv = file.path(audit_results,
                                                           "supplement_dataset_split_provenance_v7.tsv"),
  weak_dataset_paper_context_audit.tsv = file.path(audit_results,
                                                    "weak_dataset_paper_context_audit.tsv"),
  leave_one_gene_out_summary.tsv = file.path(audit_results, "leave_one_gene_out_summary.tsv"),
  gene_dominance_summary.tsv = file.path(audit_results, "gene_dominance_summary.tsv"),
  threshold_sensitivity_summary.tsv = file.path(audit_results, "threshold_sensitivity_summary.tsv"),
  leave_one_anchor_out_summary.tsv = file.path(audit_results, "leave_one_anchor_out_summary.tsv"),
  label_permutation_null_summary.tsv = file.path(extra_root, "results",
                                                  "label_permutation_null_summary.tsv"),
  baseline_signature_summary_by_group.tsv = file.path(extra_root, "results",
                                                       "baseline_signature_summary_by_group.tsv"),
  coefficient_sensitivity_summary.tsv = file.path(extra_root, "results",
                                                   "coefficient_sensitivity_summary.tsv"),
  step09_split_summary.tsv = file.path(evaluation, "step09_split_summary.tsv")
)
for (nm in names(generated_top)) {
  package_file(generated_top[[nm]], file.path(derived_root, nm), "generated",
               "Generated from the all_scored frozen-score workflow or copied original bridge analysis.")
}
package_file(file.path(ctx$paths$released_derived_root, "gene_symbol_mapping.tsv"),
             file.path(derived_root, "gene_symbol_mapping.tsv"), "released_reference",
             "Annotation mapping is a released reference input; it is not recomputed from Step09.")

generated_figure <- c(
  gene_weights.tsv = file.path(anchors, "gene_weights.tsv"),
  gene_heterogeneity.tsv = file.path(anchors, "gene_heterogeneity.tsv"),
  gene_power.tsv = file.path(anchors, "gene_power.tsv"),
  support_by_dataset.tsv = file.path(anchors, "support_by_dataset.tsv"),
  contrast_counts_by_dataset.tsv = file.path(anchors, "contrast_counts_by_dataset.tsv"),
  contrast_counts_by_dataset_6B.tsv = file.path(anchors, "contrast_counts_by_dataset_6B.tsv"),
  manuscript_dataset_role_table.tsv = file.path(audit_results, "manuscript_dataset_role_table.tsv"),
  weak_dataset_paper_context_audit.tsv = file.path(audit_results, "weak_dataset_paper_context_audit.tsv"),
  leave_one_gene_out_summary.tsv = file.path(audit_results, "leave_one_gene_out_summary.tsv"),
  gene_dominance_summary.tsv = file.path(audit_results, "gene_dominance_summary.tsv"),
  threshold_sensitivity_summary.tsv = file.path(audit_results, "threshold_sensitivity_summary.tsv"),
  threshold_sensitivity_contrast_deltas.tsv = file.path(audit_results,
                                                        "threshold_sensitivity_contrast_deltas.tsv"),
  step09_split_eval.tsv = file.path(evaluation, "step09_split_eval.tsv"),
  step09_split_sample_level.tsv = file.path(evaluation, "step09_split_sample_level.tsv"),
  label_permutation_null_summary.tsv = file.path(extra_root, "results",
                                                  "label_permutation_null_summary.tsv"),
  baseline_signature_contrast_long.tsv = file.path(extra_root, "results",
                                                    "baseline_signature_contrast_long.tsv"),
  baseline_signature_paired_contrast_comparison.tsv = file.path(extra_root, "results",
                                                                  "baseline_signature_paired_contrast_comparison.tsv"),
  baseline_signature_summary_by_group.tsv = file.path(extra_root, "results",
                                                       "baseline_signature_summary_by_group.tsv"),
  baseline_signature_scores_sample_level.tsv = file.path(extra_root, "results",
                                                          "baseline_signature_scores_sample_level.tsv"),
  coefficient_sensitivity_summary.tsv = file.path(extra_root, "results",
                                                   "coefficient_sensitivity_summary.tsv")
)
for (nm in names(generated_figure)) {
  package_file(generated_figure[[nm]], file.path(figure_root, nm), "generated",
               "Regenerated Layer2 figure input.")
}
for (nm in c("gene_symbol_mapping.tsv", "dataset_classification_authoritative.tsv",
             "five_anchor_leave_one_anchor_out_summary.tsv",
             "five_anchor_leave_one_anchor_out_contrast_details.tsv")) {
  package_file(file.path(ctx$paths$released_derived_root, "figure_inputs", nm),
               file.path(figure_root, nm), "released_reference",
               "Required by the active Layer2 runner; no matching direct bridge generator was ported.")
}

role_path <- file.path(figure_root, "manuscript_dataset_role_table.tsv")
roles <- readr::read_tsv(role_path, show_col_types = FALSE, progress = FALSE)
required_role_labels <- c("Locked anchor", "Primary acute validation",
                          "Extended validation", "Secondary support")
missing_roles <- setdiff(required_role_labels, unique(roles$role_display))
if (length(missing_roles) > 0) {
  stop("Generated role table is missing required Layer2 role display group(s): ",
       paste(missing_roles, collapse = ", "), call. = FALSE)
}
required_primary <- c("GSE119119", "GSE139529", "GSE279743")
available_primary <- intersect(required_primary, imrs_available_datasets(ctx))
generated_primary <- unique(roles$dataset_id[roles$role_display == "Primary acute validation"])
missing_primary <- setdiff(available_primary, generated_primary)
if (length(missing_primary) > 0) {
  stop("Generated role table lacks available primary acute validation dataset(s): ",
       paste(missing_primary, collapse = ", "), call. = FALSE)
}

manifest <- dplyr::bind_rows(manifest_rows)
readr::write_tsv(manifest, manifest_path)
layer2_config <- list(
  project_root = ctx$repo_root,
  derived_input_dir = derived_root,
  derived_data_dir = derived_root,
  figure_input_dir = figure_root,
  results_dir = file.path(bridge_root, "manuscript_outputs"),
  logs_dir = file.path(bridge_root, "manuscript_outputs", "logs"),
  manifests_dir = file.path(bridge_root, "manuscript_outputs", "manifests"),
  figures_dir = file.path(bridge_root, "manuscript_outputs", "figures"),
  supplementary_figures_dir = file.path(bridge_root, "manuscript_outputs", "supplementary_figures"),
  supplementary_tables_dir = file.path(bridge_root, "manuscript_outputs", "supplementary_tables"),
  manuscript_output_dir = file.path(bridge_root, "manuscript_outputs"),
  priority3_enrichment_dir = file.path(bridge_root, "manuscript_outputs",
                                       "priority3_gene_program_enrichment"),
  execute_active_scripts = FALSE,
  run_internal_readiness = FALSE,
  figure_helper_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                   "figure_helpers_v6.R"),
  panel_builder_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                   "panel_builders_v6.R"),
  workflow_diagram_script = file.path(ctx$repo_root, "scripts", "active_manuscript", "lib",
                                      "merged_workflow_v6.R"),
  frozen_gene_weights = file.path(derived_root, "frozen_gene_weights.tsv"),
  gene_symbol_mapping = file.path(derived_root, "gene_symbol_mapping.tsv"),
  gene_power = file.path(derived_root, "gene_power.tsv"),
  provenance_table = file.path(derived_root, "supplement_dataset_split_provenance_v7.tsv"),
  manuscript_role_table = file.path(figure_root, "manuscript_dataset_role_table.tsv"),
  boundary_audit_table = file.path(derived_root, "weak_dataset_paper_context_audit.tsv"),
  label_permutation_summary = file.path(derived_root, "label_permutation_null_summary.tsv"),
  leave_one_gene_out_summary = file.path(derived_root, "leave_one_gene_out_summary.tsv"),
  gene_dominance_summary = file.path(derived_root, "gene_dominance_summary.tsv"),
  threshold_sensitivity_summary = file.path(derived_root, "threshold_sensitivity_summary.tsv"),
  leave_one_anchor_out_summary = file.path(derived_root, "leave_one_anchor_out_summary.tsv"),
  comparator_benchmarking_summary = file.path(derived_root, "baseline_signature_summary_by_group.tsv"),
  coefficient_sensitivity_summary = file.path(derived_root, "coefficient_sensitivity_summary.tsv"),
  step09_split_eval = file.path(figure_root, "step09_split_eval.tsv"),
  step09_split_summary = file.path(derived_root, "step09_split_summary.tsv"),
  source_v5_dir = derived_root
)
yaml::write_yaml(layer2_config, layer2_config_path)
log_msg("Layer2 generated input package completed: ", bridge_root)
log_msg("Released data/derived tables were read as references only and were not overwritten.")
log_msg("Layer2 dry-run config: ", layer2_config_path)
