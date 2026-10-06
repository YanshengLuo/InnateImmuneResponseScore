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

s2_stem <- file.path(
  repo_root, "results_release_templates", "figures",
  "FigureS2_gene_program_enrichment_combined"
)
s2_paths <- paste0(s2_stem, c(".png", ".pdf", ".svg"))
released_s2_stem <- file.path(
  repo_root, "publication_outputs", "figures",
  "FigureS2_gene_program_enrichment_combined"
)
for (extension in c(".png", ".pdf", ".svg")) {
  target <- paste0(s2_stem, extension)
  released <- paste0(released_s2_stem, extension)
  if (!file.exists(target) || file.info(target)$size <= 0) {
    if (!file.exists(released) || file.info(released)$size <= 0) {
      stop("Missing released Supplementary Figure S2 snapshot: ", released,
           call. = FALSE)
    }
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(released, target, overwrite = TRUE)) {
      stop("Could not restore released Supplementary Figure S2: ", target,
           call. = FALSE)
    }
  }
}
missing_s2 <- s2_paths[!file.exists(s2_paths) | file.info(s2_paths)$size <= 0]
if (length(missing_s2) > 0L) {
  stop("Missing or empty retained Supplementary Figure S2 output(s): ",
       paste(missing_s2, collapse = "; "), call. = FALSE)
}
s2_svg_text <- paste(readLines(paste0(s2_stem, ".svg"), warn = FALSE,
                               encoding = "UTF-8"), collapse = "\n")
s2_required <- c("GO Biological Process", "Reactome", "MSigDB Hallmark")
s2_missing_text <- s2_required[!vapply(
  s2_required, grepl, logical(1), x = s2_svg_text, fixed = TRUE
)]
if (length(s2_missing_text) > 0L) {
  stop("Retained Supplementary Figure S2 is missing required program labels: ",
       paste(s2_missing_text, collapse = "; "), call. = FALSE)
}
message("Verified retained Supplementary Figure S2 without rerunning enrichment: ",
        paste(basename(s2_paths), collapse = ", "))

default_rplots <- file.path(repo_root, "Rplots.pdf")
if (file.exists(default_rplots)) {
  grDevices::graphics.off()
  unlink(default_rplots)
  if (file.exists(default_rplots)) {
    stop("Could not remove unintended default graphics device output: ",
         default_rplots, call. = FALSE)
  }
  message("Removed unintended default graphics device output: ", default_rplots)
}
