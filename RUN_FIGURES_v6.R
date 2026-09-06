#!/usr/bin/env Rscript

# One-click manuscript figure regeneration for the extracted IMRS reviewer package.
# Safe to run from RStudio (Source or Run) or with Rscript; no setwd() is required.

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

this_file <- imrs_detect_script_path("RUN_FIGURES_v6.R")
starts <- unique(c(if (!is.na(this_file)) dirname(this_file) else character(), getwd()))
repo_root <- NA_character_
for (start in starts) {
  candidate <- imrs_find_repo_root_bootstrap(start)
  if (!is.na(candidate)) { repo_root <- candidate; break }
}
if (is.na(repo_root)) {
  stop("Could not locate the extracted IMRS repository root.", call. = FALSE)
}

figure_script <- file.path(repo_root, "scripts", "active_manuscript", "00_generate_manuscript_figures_v6.R")
if (!file.exists(figure_script)) stop("Missing figure entry script: ", figure_script, call. = FALSE)

Sys.setenv(IMRS_REPOSITORY_ROOT = repo_root)
message("IMRS repository root: ", repo_root)
message("Running manuscript figure generation: ", figure_script)
source(figure_script, chdir = FALSE)
