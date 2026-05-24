#!/usr/bin/env Rscript

# Reviewer-facing manuscript-output runner for IMRS v6.
# This runner executes only scripts in scripts/active_manuscript.
# It intentionally excludes raw-data retrieval, HiPerGator/HPC scripts,
# archive_legacy scripts, and full-pipeline reconstruction steps.
# Active manuscript reproduction starts from released frozen derived inputs,
# including data_release_templates/derived/frozen_gene_weights.tsv and
# data_release_templates/derived/gene_power.tsv. Full frozen-gene reconstruction
# from clean-gene / anchor inputs is documented separately under
# scripts/full_pipeline/frozen_gene_reconstruction and is not run here.

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

read_simple_yaml <- function(file) {
  if (!file.exists(file)) stop("Config file does not exist: ", file, call. = FALSE)
  lines <- readLines(file, warn = FALSE)
  lines <- lines[grepl(":", lines) & !grepl("^\\s*#", lines)]
  out <- list()
  for (line in lines) {
    key <- trimws(sub(":.*$", "", line))
    value <- trimws(sub("^[^:]+:", "", line))
    value <- trimws(gsub("^['\"]|['\"]$", "", value))
    if (tolower(value) %in% c("true", "false")) value <- tolower(value) == "true"
    out[[key]] <- value
  }
  out
}

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                      error = function(e) NA_character_)
if (is.na(this_file)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
}
repo_root <- if (!is.na(this_file)) {
  normalizePath(dirname(this_file), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
if (!dir.exists(file.path(repo_root, "config"))) {
  stop("Run this script from the IMRS_Repository_Scripts root or keep it in that root.", call. = FALSE)
}

config_file <- file.path(repo_root, "config", "config.yml")
if (!file.exists(config_file)) config_file <- file.path(repo_root, "config", "config_template.yml")
config <- read_simple_yaml(config_file)

resolve_path <- function(value, default = NULL) {
  if (is.null(value) || is.na(value) || !nzchar(value)) value <- default
  if (is.null(value) || is.na(value) || !nzchar(value)) return(NA_character_)
  out <- if (grepl("^([A-Za-z]:|/)", value)) value else file.path(repo_root, value)
  normalizePath(out, winslash = "/", mustWork = FALSE)
}

project_root <- resolve_path(config$project_root %||% ".")
logs_dir <- resolve_path(config$logs_dir %||% "results_release_templates/logs")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs_dir, paste0("run_all_manuscript_outputs_v6_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
  invisible(msg)
}

path_exists <- function(value, default = NULL) {
  p <- resolve_path(value, default)
  !is.na(p) && file.exists(p)
}

script_path <- function(name) file.path(repo_root, "scripts", "active_manuscript", name)
required_path <- function(field, default = NULL) resolve_path(config[[field]], default)

steps <- data.frame(
  order = seq_len(5),
  step = c(
    "main_figure_regeneration",
    "targeted_v6_figure_regeneration",
    "supplementary_figure_s2_enrichment",
    "supplementary_tables_s1_s5",
    "nar_gb_readiness"
  ),
  script = c(
    script_path("00_generate_all_reorganized_figures_revised_v5.R"),
    script_path("01_regenerate_changed_v6.R"),
    script_path("run_IMRS_retained_gene_enrichment_v6.R"),
    script_path("build_supplementary_tables_v6.R"),
    script_path("build_NAR_GB_readiness_v6.R")
  ),
  expected_output = c(
    "Main Figures 1-5 and Supplementary Figure S1 outputs in manuscript_output_dir.",
    "Targeted v6 Figure 1 and revised panel outputs in manuscript_output_dir.",
    "Supplementary Figure S2 and retained-gene enrichment tables in priority3_enrichment_dir.",
    "Supplementary Tables S1-S5 in supplementary_tables_dir.",
    "NAR G&B readiness files in nar_readiness_dir."
  ),
  required_or_optional = c("optional", "optional", "required", "required", "optional"),
  stringsAsFactors = FALSE
)

input_checks <- list(
  main_figure_regeneration = c("legacy_revised_plots_v2_dir", "legacy_revised_plots_v3_dir", "legacy_revised_plots_v4_dir", "figure_helper_script"),
  targeted_v6_figure_regeneration = c("legacy_figure_generator_script", "figure_helper_script", "manuscript_output_dir"),
  supplementary_figure_s2_enrichment = c("frozen_gene_weights", "gene_power"),
  supplementary_tables_s1_s5 = c("provenance_table", "boundary_audit_table", "label_permutation_summary", "leave_one_gene_out_summary", "gene_dominance_summary", "threshold_sensitivity_summary", "leave_one_anchor_out_summary"),
  nar_gb_readiness = c("supplementary_tables_dir", "manuscript_output_dir")
)

check_step_inputs <- function(step_name) {
  fields <- input_checks[[step_name]]
  if (is.null(fields)) return(data.frame(field = character(), path = character(), exists = logical()))
  data.frame(
    field = fields,
    path = vapply(fields, function(x) required_path(x), character(1)),
    exists = vapply(fields, function(x) path_exists(config[[x]]), logical(1)),
    stringsAsFactors = FALSE
  )
}

if (any(grepl("scripts/(hpc_hypergator|archive_legacy)", normalizePath(steps$script, winslash = "/", mustWork = FALSE)))) {
  stop("Runner configuration attempted to call hpc_hypergator or archive_legacy scripts.", call. = FALSE)
}

Sys.setenv(
  IMRS_REPOSITORY_ROOT = repo_root,
  IMRS_PROJECT_ROOT = project_root,
  IMRS_MANUSCRIPT_OUTPUT_DIR = required_path("manuscript_output_dir")
)

log_msg("IMRS v6 manuscript-output runner started.")
log_msg("Repository root: ", repo_root)
log_msg("Config file: ", config_file)
log_msg("Project root: ", project_root)
log_msg("execute_active_scripts: ", config$execute_active_scripts %||% FALSE)

find_rscript <- function() {
  candidates <- unique(c(
    if (.Platform$OS.type == "windows") file.path(R.home("bin"), "Rscript.exe") else file.path(R.home("bin"), "Rscript"),
    file.path(R.home("bin"), "Rscript"),
    Sys.which("Rscript")
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop("Could not locate Rscript. Checked: ", paste(candidates, collapse = "; "), call. = FALSE)
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

steps$script_exists <- file.exists(steps$script)
steps$missing_inputs <- vapply(steps$step, function(step_name) {
  checks <- check_step_inputs(step_name)
  paste(checks$field[!checks$exists], collapse = "; ")
}, character(1))
steps$can_run_now <- steps$script_exists & !nzchar(steps$missing_inputs)
steps$status <- ifelse(steps$can_run_now, "ready", ifelse(steps$required_or_optional == "optional", "optional_blocked", "required_blocked"))

write.table(steps, file = file.path(logs_dir, "run_all_manuscript_outputs_v6_steps.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

blocked_required <- steps[steps$status == "required_blocked", , drop = FALSE]
if (nrow(blocked_required) > 0) {
  msg <- paste0(blocked_required$step, " missing: ", blocked_required$missing_inputs, collapse = " | ")
  log_msg("Required active step inputs are missing: ", msg)
  if (isTRUE(config$execute_active_scripts)) {
    stop("Required active step inputs are missing: ", msg, call. = FALSE)
  }
}

if (!isTRUE(config$execute_active_scripts)) {
  log_msg("execute_active_scripts is false. Wrote step checklist only; no outputs regenerated.")
  print(steps[, c("order", "step", "required_or_optional", "status", "missing_inputs", "expected_output")])
  quit(status = 0)
}

run_one_step <- function(step_row) {
  log_msg("Starting step ", step_row$order, ": ", step_row$step)
  rscript <- find_rscript()
  out_log <- file.path(logs_dir, paste0(step_row$step, "_stdout_stderr.log"))
  dir.create(dirname(out_log), recursive = TRUE, showWarnings = FALSE)
  log_msg("Step Rscript: ", rscript)
  log_msg("Step script: ", step_row$script)
  log_msg("Step stdout/stderr log: ", out_log)
  Sys.setenv(
    IMRS_REPOSITORY_ROOT = repo_root,
    IMRS_PROJECT_ROOT = project_root,
    IMRS_MANUSCRIPT_OUTPUT_DIR = required_path("manuscript_output_dir")
  )
  status <- system2(rscript,
                    args = shQuote(step_row$script),
                    stdout = out_log,
                    stderr = out_log)
  if (is.null(status)) status <- 0L
  status <- as.integer(status)
  log_msg("Step exit status ", step_row$step, ": ", status)
  if (!identical(status, 0L)) {
    stop("Step failed: ", step_row$step, " with exit status ", status, ". See log: ", out_log, call. = FALSE)
  }
  log_msg("Completed step ", step_row$order, ": ", step_row$step)
}

executed <- character()
skipped <- character()
for (i in seq_len(nrow(steps))) {
  row <- steps[i, , drop = FALSE]
  if (!row$can_run_now) {
    if (row$required_or_optional == "optional") {
      log_msg("Skipping optional step ", row$step, "; missing inputs: ", row$missing_inputs)
      skipped <- c(skipped, row$step)
      next
    }
    stop("Required step cannot run: ", row$step, "; missing inputs: ", row$missing_inputs, call. = FALSE)
  }
  run_one_step(row)
  executed <- c(executed, row$step)
}

summary_lines <- c(
  "IMRS v6 active manuscript run summary",
  paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Executed steps: ", ifelse(length(executed) == 0, "none", paste(executed, collapse = ", "))),
  paste0("Skipped optional steps: ", ifelse(length(skipped) == 0, "none", paste(skipped, collapse = ", "))),
  paste0("Log file: ", log_file)
)
writeLines(summary_lines, file.path(logs_dir, "run_all_manuscript_outputs_v6_summary.txt"), useBytes = TRUE)
log_msg("IMRS v6 manuscript-output runner finished.")
cat(paste(summary_lines, collapse = "\n"), "\n")
