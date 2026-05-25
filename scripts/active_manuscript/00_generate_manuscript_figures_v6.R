#!/usr/bin/env Rscript

# IMRS v6 reviewer-facing figure-generation/assembly layer.
# Reads released derived figure inputs and writes only to configured repository outputs.

options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                      error = function(e) NA_character_)
if (is.na(this_file)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
}
active_config_helper <- file.path(dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE)),
                                  "lib", "active_config.R")
if (!file.exists(active_config_helper)) {
  active_config_helper <- file.path(getwd(), "scripts", "active_manuscript", "lib", "active_config.R")
}
source(active_config_helper)

config <- imrs_load_active_config(dirname(active_config_helper))
project_root <- imrs_project_root(config)
input_root <- imrs_config_field_path(config, "figure_input_dir")
v5_root <- imrs_config_field_path(config, "figures_dir")
helper <- imrs_config_field_path(config, "figure_helper_script",
                                 "scripts/active_manuscript/lib/figure_helpers_v6.R")
panel_builder <- imrs_config_field_path(config, "panel_builder_script",
                                        "scripts/active_manuscript/lib/panel_builders_v6.R")
workflow_script <- imrs_config_field_path(config, "workflow_diagram_script",
                                          "scripts/active_manuscript/lib/merged_workflow_v6.R")

if (!dir.exists(input_root)) {
  stop("Released derived figure-input directory does not exist: ", input_root, call. = FALSE)
}
required_scripts <- c(helper, panel_builder, workflow_script)
missing_scripts <- required_scripts[!file.exists(required_scripts)]
if (length(missing_scripts) > 0) {
  stop("Missing repo-contained figure implementation script(s): ",
       paste(missing_scripts, collapse = "; "), call. = FALSE)
}

source(helper)
stop_if_missing_packages_v5()

dir.create(v5_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(v5_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(v5_root, "intermediate_panels"), recursive = TRUE, showWarnings = FALSE)

baseline_inputs <- newest_file_snapshot_v5(input_root)

log_msg_v5("Starting IMRS v6 generation.")
log_msg_v5("Released figure-input root: ", input_root)
log_msg_v5("Output root: ", v5_root)
Sys.setenv(IMRS_FIGURE_INPUT_DIR = input_root, IMRS_MANUSCRIPT_OUTPUT_DIR = v5_root)

panel_manifest <- render_v5_panels(
  project_root = project_root,
  input_root = input_root,
  v5_root = v5_root,
  dpi = 400,
  panel_builder_script = panel_builder,
  workflow_script = workflow_script
)

figure_manifest <- assemble_v5_figures(
  v5_root = v5_root,
  panel_manifest = panel_manifest,
  dpi = 400
)

manifest <- manifest_long_v5(panel_manifest, figure_manifest, v5_root)
wording_audit <- wording_audit_rows_v5()
write_tsv_v5(wording_audit, file.path(v5_root, "tables", "v5_wording_audit.tsv"))
checks <- validate_v5_outputs(v5_root, manifest, baseline_inputs, baseline_inputs, baseline_inputs)
log_path <- write_generation_log_v5(v5_root, manifest$file_path, checks)

cat("\nIMRS v6 generation complete\n")
cat("---------------------------\n")
cat("Total PNG files generated: ", sum(manifest$format == "png"), "\n", sep = "")
cat("Total PDF files generated: ", sum(manifest$format == "pdf"), "\n", sep = "")
cat("Total SVG files generated: ", sum(manifest$format == "svg"), "\n", sep = "")
cat("Output folder: ", v5_root, "\n", sep = "")
cat("Skipped panels: ", ifelse(length(v5_skipped) == 0, "none", paste(v5_skipped, collapse = ", ")), "\n", sep = "")
cat("Warnings: ", ifelse(length(v5_warnings) == 0, "none", paste(v5_warnings, collapse = " | ")), "\n", sep = "")
cat("Confirmation: released derived input tables were read without modification; outputs were written to the configured repository results folder.\n")
cat("Manifest: ", file.path(v5_root, "figure_v5_manifest.tsv"), "\n", sep = "")
cat("Wording audit: ", file.path(v5_root, "tables", "v5_wording_audit.tsv"), "\n", sep = "")
cat("Log: ", log_path, "\n", sep = "")
