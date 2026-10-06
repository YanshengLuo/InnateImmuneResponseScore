# ---------------------------------------------------------------------------
# manuscript_universe.R
#
# THE SINGLE AUTHORITATIVE DEFINITION of which datasets contribute to scored
# manuscript outputs.
#
# Every scored manuscript-facing table and figure must pass its rows through
# imrs_filter_scored() exactly once. No other file may gate a scored output on
# a literal dataset identifier.
#
# The excluded arms are deliberately NOT deleted from metadata, provenance or
# split-design artifacts. They remain fully documented; they are only barred
# from scored aggregates.
#
# Source of truth: data/curated_metadata/manuscript_scored_universe.tsv
# ---------------------------------------------------------------------------

imrs_universe_path <- function(project_root = ".") {
  file.path(project_root, "data", "curated_metadata", "manuscript_scored_universe.tsv")
}

imrs_load_universe <- function(project_root = ".") {
  p <- imrs_universe_path(project_root)
  if (!file.exists(p)) {
    stop("manuscript_scored_universe.tsv not found at: ", p,
         "\nThis file defines the scored manuscript universe and is required.",
         call. = FALSE)
  }
  u <- readr::read_tsv(p, col_types = readr::cols(.default = readr::col_character()),
                       progress = FALSE)
  required <- c("dataset_id", "in_scored_manuscript_universe",
                "retained_in_metadata_and_provenance", "exclusion_reason")
  missing <- setdiff(required, names(u))
  if (length(missing) > 0) {
    stop("manuscript_scored_universe.tsv is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  u$in_scored <- toupper(trimws(u$in_scored_manuscript_universe)) == "TRUE"
  u$retained <- toupper(trimws(u$retained_in_metadata_and_provenance)) == "TRUE"
  if (!all(u$retained)) {
    stop("Every dataset must be retained in metadata/provenance; offending id(s): ",
         paste(u$dataset_id[!u$retained], collapse = ", "), call. = FALSE)
  }
  u
}

#' Dataset ids that contribute to scored manuscript outputs.
imrs_scored_ids <- function(project_root = ".") {
  u <- imrs_load_universe(project_root)
  u$dataset_id[u$in_scored]
}

#' Dataset ids retained in metadata but excluded from scored outputs.
imrs_excluded_ids <- function(project_root = ".") {
  u <- imrs_load_universe(project_root)
  u$dataset_id[!u$in_scored]
}

#' Filter a table to the scored manuscript universe.
#'
#' Fails loudly on any dataset id that the manifest does not declare, so a new
#' dataset can never slip silently into (or out of) a scored output.
#'
#' @param df       data frame carrying a dataset identifier column
#' @param id_col   name of that column (default "dataset_id")
#' @param label    short description used in the console message
imrs_filter_scored <- function(df, id_col = "dataset_id", label = "table",
                               project_root = ".") {
  if (!(id_col %in% names(df))) {
    stop("imrs_filter_scored(): column '", id_col, "' not present in ", label,
         ". Available: ", paste(names(df), collapse = ", "), call. = FALSE)
  }
  u <- imrs_load_universe(project_root)
  ids <- unique(trimws(as.character(df[[id_col]])))
  undeclared <- setdiff(ids[nzchar(ids) & !is.na(ids)], u$dataset_id)
  if (length(undeclared) > 0) {
    stop("imrs_filter_scored(): dataset id(s) not declared in ",
         "manuscript_scored_universe.tsv: ", paste(undeclared, collapse = ", "),
         "\nAdd them to the manifest before they may appear in a scored output.",
         call. = FALSE)
  }
  keep_ids <- u$dataset_id[u$in_scored]
  before <- nrow(df)
  out <- df[trimws(as.character(df[[id_col]])) %in% keep_ids, , drop = FALSE]
  dropped <- before - nrow(out)
  if (dropped > 0) {
    message(sprintf(
      "[manuscript universe] %s: %d -> %d rows (%d excluded-arm row(s) removed; retained in metadata)",
      label, before, nrow(out), dropped))
  }
  out
}

#' De-duplicate split-level provenance.
#'
#' The GSE262515 splits are registered twice in the split-design tree: once
#' under the generic GSE262515_design/ directory and again under the
#' arm-specific GSE262515_cell_line_design/ and GSE262515_tissue_design/
#' directories. Both registrations normalise to the same (dataset_id, split_id)
#' with identical payloads, which inflated per-dataset counts in Supplementary
#' Table S1 (GSE262515_tissue reported 6 positive contrasts out of 3).
imrs_dedup_provenance <- function(df, label = "provenance") {
  if (!all(c("dataset_id", "split_id") %in% names(df))) return(df)
  before <- nrow(df)
  out <- df[!duplicated(df[, c("dataset_id", "split_id")]), , drop = FALSE]
  if (before != nrow(out)) {
    message(sprintf("[manuscript universe] %s de-duplicated: %d -> %d rows",
                    label, before, nrow(out)))
  }
  out
}
