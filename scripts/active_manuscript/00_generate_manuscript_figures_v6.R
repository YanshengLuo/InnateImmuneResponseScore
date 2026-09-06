#!/usr/bin/env Rscript

# IMRS v6 reviewer-facing figure-generation/assembly layer.
# Reads released derived figure inputs and writes only to configured repository outputs.

options(stringsAsFactors = FALSE)

imrs_detect_script_path <- function(expected_basename = NULL) {
  candidates <- character()

  # Rscript --file=... invocation.
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    candidates <- c(candidates, sub("^--file=", "", file_arg[[1L]]))
  }

  # source()/sys.source() invocation (search all active frames, not only frame 1).
  frame_paths <- unlist(lapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.null(path) || length(path) == 0L) character() else as.character(path[[1L]])
  }), use.names = FALSE)
  candidates <- c(candidates, frame_paths)

  # RStudio "Run" / "Source" invocation, including running selected lines.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    editor_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    candidates <- c(candidates, editor_path)
  }

  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (length(candidates) > 0L) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
    if (!is.null(expected_basename) && nzchar(expected_basename)) {
      matching <- candidates[basename(candidates) == expected_basename]
      if (length(matching) > 0L) candidates <- c(matching, candidates)
    }
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
  }

  # Last-resort support when the working directory is the script directory.
  if (!is.null(expected_basename) && nzchar(expected_basename)) {
    cwd_candidate <- file.path(getwd(), expected_basename)
    if (file.exists(cwd_candidate)) {
      return(normalizePath(cwd_candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

imrs_find_repo_root_bootstrap <- function(start = getwd()) {
  if (is.null(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
    start <- getwd()
  }
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)
  if (file.exists(current) && !dir.exists(current)) current <- dirname(current)

  repeat {
    marker_config <- file.path(current, "config", "config_template.yml")
    marker_active <- file.path(current, "scripts", "active_manuscript", "lib", "active_config.R")
    if (file.exists(marker_config) && file.exists(marker_active)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  NA_character_
}

this_file <- imrs_detect_script_path("00_generate_manuscript_figures_v6.R")
bootstrap_starts <- unique(c(
  if (!is.na(this_file)) dirname(this_file) else character(),
  getwd()
))
repo_root_bootstrap <- NA_character_
for (start in bootstrap_starts) {
  candidate_root <- imrs_find_repo_root_bootstrap(start)
  if (!is.na(candidate_root)) {
    repo_root_bootstrap <- candidate_root
    break
  }
}
if (is.na(repo_root_bootstrap)) {
  stop(
    "Could not locate the IMRS repository root. Keep the extracted repository structure intact and run this file with RStudio Source/Run or Rscript.",
    call. = FALSE
  )
}
active_config_helper <- file.path(
  repo_root_bootstrap, "scripts", "active_manuscript", "lib", "active_config.R"
)
source(active_config_helper)


config <- imrs_load_active_config(repo_root_bootstrap)
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
