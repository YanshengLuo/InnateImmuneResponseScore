#!/usr/bin/env Rscript

# Reviewer-facing Layer 2 manuscript-output runner for IMRS v6.
# It starts from released derived inputs and runs only repo-contained active
# manuscript scripts. It does not run raw-data retrieval, HPC jobs, full
# framework reconstruction, or frozen-gene reconstruction.

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a) || !nzchar(as.character(a))) b else a

read_simple_yaml <- function(file) {
  if (!file.exists(file)) stop("Config file does not exist: ", file, call. = FALSE)
  lines <- readLines(file, warn = FALSE)
  lines <- lines[grepl(":", lines) & !grepl("^\\s*#", lines)]
  out <- list()
  for (line in lines) {
    key <- trimws(sub(":.*$", "", line))
    value <- trimws(sub("^[^:]+:", "", line))
    value <- trimws(gsub("^['\"]|['\"]$", "", value))
    if (tolower(value) %in% c("true", "false", "yes", "no")) {
      value <- tolower(value) %in% c("true", "yes")
    }
    out[[key]] <- value
  }
  out
}

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

runner_file <- imrs_detect_script_path("run_all_manuscript_outputs_v6.R")
bootstrap_starts <- unique(c(
  if (!is.na(runner_file)) dirname(runner_file) else character(),
  getwd()
))
repo_root <- NA_character_
for (start in bootstrap_starts) {
  candidate_root <- imrs_find_repo_root_bootstrap(start)
  if (!is.na(candidate_root)) {
    repo_root <- candidate_root
    break
  }
}
if (is.na(repo_root)) {
  stop(
    "Could not locate the IMRS repository root. Keep the extracted repository structure intact and run this file with RStudio Source/Run or Rscript.",
    call. = FALSE
  )
}


runner_args <- commandArgs(trailingOnly = TRUE)
config_flag <- match("--config", runner_args)
if (!is.na(config_flag)) {
  if (config_flag == length(runner_args)) stop("Missing value after --config.", call. = FALSE)
  config_arg <- runner_args[[config_flag + 1L]]
  config_file <- if (grepl("^([A-Za-z]:|/)", config_arg)) config_arg else
    file.path(repo_root, config_arg)
  config_file <- normalizePath(config_file, winslash = "/", mustWork = FALSE)
} else {
  config_file <- file.path(repo_root, "config", "config.yml")
  if (!file.exists(config_file)) config_file <- file.path(repo_root, "config", "config_template.yml")
}
config <- read_simple_yaml(config_file)

resolve_path <- function(value, default = NULL) {
  value <- value %||% default
  if (is.null(value) || is.na(value) || !nzchar(value)) return(NA_character_)
  path <- if (grepl("^([A-Za-z]:|/)", value)) value else file.path(repo_root, value)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

logs_dir <- resolve_path(config$logs_dir, "results_release_templates/logs")
manifests_dir <- resolve_path(config$manifests_dir, "results_release_templates/manifests")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifests_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs_dir, paste0("run_all_manuscript_outputs_v6_",
                                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(line)
}

active_script <- function(name) file.path(repo_root, "scripts", "active_manuscript", name)
internal_script <- function(name) file.path(repo_root, "scripts", "optional_internal", name)

figure_inputs <- file.path(resolve_path(config$figure_input_dir, "data/derived/figure_inputs"), c(
  "manuscript_dataset_role_table.tsv", "dataset_classification_authoritative.tsv",
  "support_by_dataset.tsv", "contrast_counts_by_dataset.tsv", "contrast_counts_by_dataset_6B.tsv",
  "gene_weights.tsv", "gene_symbol_mapping.tsv", "gene_heterogeneity.tsv", "gene_power.tsv",
  "leave_one_gene_out_summary.tsv", "gene_dominance_summary.tsv",
  "threshold_sensitivity_summary.tsv", "threshold_sensitivity_contrast_deltas.tsv",
  "five_anchor_leave_one_anchor_out_summary.tsv",
  "five_anchor_leave_one_anchor_out_contrast_details.tsv",
  "weak_dataset_paper_context_audit.tsv", "step09_split_eval.tsv",
  "step09_split_sample_level.tsv", "label_permutation_null_summary.tsv",
  "baseline_signature_contrast_long.tsv", "baseline_signature_paired_contrast_comparison.tsv",
  "baseline_signature_summary_by_group.tsv", "baseline_signature_scores_sample_level.tsv",
  "coefficient_sensitivity_summary.tsv"
))

input_specs <- list(
  manuscript_figures = c(
    resolve_path(config$figure_helper_script, "scripts/active_manuscript/lib/figure_helpers_v6.R"),
    resolve_path(config$panel_builder_script, "scripts/active_manuscript/lib/panel_builders_v6.R"),
    resolve_path(config$workflow_diagram_script, "scripts/active_manuscript/lib/merged_workflow_v6.R"),
    figure_inputs
  ),
  priority3_gene_program_enrichment = c(
    resolve_path(config$frozen_gene_weights, "data/derived/frozen_gene_weights.tsv"),
    resolve_path(config$gene_power, "data/derived/gene_power.tsv"),
    resolve_path(config$gene_symbol_mapping, "data/derived/gene_symbol_mapping.tsv")
  ),
  supplementary_tables_s1_s5 = c(
    resolve_path(config$provenance_table, "data/derived/supplement_dataset_split_provenance_v7.tsv"),
    resolve_path(config$boundary_audit_table, "data/derived/weak_dataset_paper_context_audit.tsv"),
    resolve_path(config$label_permutation_summary, "data/derived/label_permutation_null_summary.tsv"),
    resolve_path(config$leave_one_gene_out_summary, "data/derived/leave_one_gene_out_summary.tsv"),
    resolve_path(config$gene_dominance_summary, "data/derived/gene_dominance_summary.tsv"),
    resolve_path(config$threshold_sensitivity_summary, "data/derived/threshold_sensitivity_summary.tsv"),
    resolve_path(config$leave_one_anchor_out_summary, "data/derived/leave_one_anchor_out_summary.tsv"),
    resolve_path(config$comparator_benchmarking_summary, "data/derived/baseline_signature_summary_by_group.tsv"),
    resolve_path(config$coefficient_sensitivity_summary, "data/derived/coefficient_sensitivity_summary.tsv")
  ),
  internal_readiness = c(
    resolve_path(config$supplementary_tables_dir, "results_release_templates/supplementary_tables")
  )
)

run_internal <- isTRUE(config$run_internal_readiness)
steps <- data.frame(
  order = seq_len(4),
  step = c("manuscript_figures", "priority3_gene_program_enrichment",
           "supplementary_tables_s1_s5", "internal_readiness"),
  script = c(active_script("00_generate_manuscript_figures_v6.R"),
             active_script("02_run_priority3_gene_program_enrichment_v6.R"),
             active_script("01_build_supplementary_tables_v6.R"),
             internal_script("build_NAR_GB_readiness_v6.R")),
  expected_output = c(
    "Figures 1-5 and supplementary manuscript figure assemblies in results_release_templates/figures.",
    "Supplementary Figure S2 and retained-gene enrichment tables in results_release_templates/priority3_gene_program_enrichment.",
    "Supplementary Tables S1-S5 in results_release_templates/supplementary_tables.",
    "Optional internal readiness documentation in results_release_templates/internal_qc/NAR_GB_readiness."
  ),
  required_or_optional = c("required", "required", "required", "optional"),
  enabled = c(TRUE, TRUE, TRUE, run_internal),
  stringsAsFactors = FALSE
)

forbidden <- grepl("scripts/(hpc_hypergator|archive_legacy|full_pipeline)",
                   normalizePath(steps$script, winslash = "/", mustWork = FALSE))
if (any(forbidden)) {
  stop("Runner configuration attempted to call a forbidden pipeline/HPC/archive script: ",
       paste(steps$script[forbidden], collapse = "; "), call. = FALSE)
}

steps$script_exists <- file.exists(steps$script)
steps$missing_inputs <- vapply(steps$step, function(step_name) {
  paths <- input_specs[[step_name]]
  paste(paths[!file.exists(paths)], collapse = "; ")
}, character(1))
steps$can_run_now <- steps$enabled & steps$script_exists & !nzchar(steps$missing_inputs)
steps$status <- ifelse(!steps$enabled, "disabled_by_config",
                       ifelse(steps$can_run_now, "ready",
                              ifelse(steps$required_or_optional == "required",
                                     "required_blocked", "optional_blocked")))

checklist_path <- file.path(manifests_dir, "active_manuscript_step_checklist.tsv")
write.table(steps, checklist_path, sep = "\t", row.names = FALSE, quote = FALSE)

Sys.setenv(IMRS_REPOSITORY_ROOT = repo_root,
           IMRS_PROJECT_ROOT = repo_root,
           IMRS_FIGURE_INPUT_DIR = resolve_path(config$figure_input_dir, "data/derived/figure_inputs"),
           IMRS_ACTIVE_CONFIG = normalizePath(config_file, winslash = "/", mustWork = FALSE))

log_msg("IMRS Layer 2 reviewer-facing manuscript-output runner started.")
log_msg("Repository root: ", repo_root)
log_msg("Config file: ", config_file)
log_msg("execute_active_scripts: ", config$execute_active_scripts %||% FALSE)
log_msg("run_internal_readiness: ", run_internal)
log_msg("Checklist: ", checklist_path)

blocked_required <- steps[steps$status == "required_blocked", , drop = FALSE]
if (nrow(blocked_required) > 0) {
  message_text <- paste0(blocked_required$step, " missing: ", blocked_required$missing_inputs,
                         collapse = " | ")
  log_msg("Required Layer 2 step blocked: ", message_text)
  if (isTRUE(config$execute_active_scripts)) stop(message_text, call. = FALSE)
}

if (!isTRUE(config$execute_active_scripts)) {
  log_msg("execute_active_scripts is false. Complete checklist written; no active script executed.")
  print(steps[, c("order", "step", "required_or_optional", "status",
                  "missing_inputs", "expected_output")])
  quit(status = 0)
}

find_rscript <- function() {
  candidates <- unique(c(
    if (.Platform$OS.type == "windows") file.path(R.home("bin"), "Rscript.exe") else file.path(R.home("bin"), "Rscript"),
    Sys.which("Rscript")
  ))
  hits <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(hits)) stop("Could not locate Rscript executable.", call. = FALSE)
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

run_one_step <- function(row) {
  rscript <- find_rscript()
  step_log <- file.path(logs_dir, paste0(row$step, "_stdout_stderr.log"))
  log_msg("Executing step ", row$step, " using ", row$script)
  status <- system2(rscript, args = shQuote(row$script), stdout = step_log, stderr = step_log)
  status <- if (is.null(status)) 0L else as.integer(status)
  log_msg("Step ", row$step, " exit status: ", status, "; log: ", step_log)
  if (status != 0L) stop("Step failed: ", row$step, ". See ", step_log, call. = FALSE)
}

executed <- character()
skipped <- character()
for (i in seq_len(nrow(steps))) {
  row <- steps[i, , drop = FALSE]
  if (!row$enabled || !row$can_run_now) {
    if (row$required_or_optional == "optional") {
      skipped <- c(skipped, row$step)
      log_msg("Skipping optional step ", row$step, " (status: ", row$status, ").")
      next
    }
    stop("Required step cannot run: ", row$step, "; missing inputs: ", row$missing_inputs, call. = FALSE)
  }
  run_one_step(row)
  executed <- c(executed, row$step)
}

summary_lines <- c(
  "IMRS Layer 2 reviewer-facing manuscript-output run summary",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Executed steps: ", paste(executed, collapse = ", ")),
  paste0("Skipped optional steps: ", ifelse(length(skipped), paste(skipped, collapse = ", "), "none")),
  "Excluded from this run: raw-data retrieval, HPC/SLURM, full pipeline, frozen-gene reconstruction.",
  paste0("Checklist: ", checklist_path),
  paste0("Log file: ", log_file)
)
summary_path <- file.path(logs_dir, "run_all_manuscript_outputs_v6_summary.txt")
writeLines(summary_lines, summary_path, useBytes = TRUE)
log_msg("IMRS Layer 2 manuscript-output runner finished.")
cat(paste(summary_lines, collapse = "\n"), "\n")
