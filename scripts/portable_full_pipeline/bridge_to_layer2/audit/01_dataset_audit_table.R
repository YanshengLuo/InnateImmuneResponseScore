#!/usr/bin/env Rscript

imrs_detect_script_path_local <- function(expected_basename = NULL) {
  candidates <- character()
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) candidates <- c(candidates, sub("^--file=", "", file_arg[[1L]]))
  frame_paths <- unlist(lapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.null(path) || length(path) == 0L) character() else as.character(path[[1L]])
  }), use.names = FALSE)
  candidates <- c(candidates, frame_paths)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    candidates <- c(candidates, tryCatch(rstudioapi::getActiveDocumentContext()$path,
                                         error = function(e) ""))
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
  if (!is.null(expected_basename) && nzchar(expected_basename)) {
    candidate <- file.path(getwd(), expected_basename)
    if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

this_script <- imrs_detect_script_path_local("01_dataset_audit_table.R")
if (is.na(this_script)) stop("Could not identify this script path. Use RStudio Source/Run or Rscript.", call. = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "01_dataset_audit_table"

get_value_from_split <- function(split, condition, candidates) {
  if (nrow(split) == 0) return(NA_character_)
  col <- first_present_col(split, candidates)
  if (is.na(col)) return(NA_character_)
  vals <- split %>%
    filter(condition_simple == condition) %>%
    pull(.data[[col]])
  collapse_unique_chr(vals)
}

get_any_value <- function(split, candidates) {
  col <- first_present_col(split, candidates)
  if (is.na(col)) return(NA_character_)
  collapse_unique_chr(split[[col]])
}

make_inclusion_rationale <- function(role, time_group) {
  case_when(
    role == "anchor" ~ "Strict anchor used for portability checks and frozen-score validation.",
    role == "calibration" & time_group == "acute_main" ~
      "Non-strict acute calibration dataset; supports external acute validation without expanding anchor set.",
    role == "calibration" & time_group == "extended_transfer" ~
      "Calibration dataset beyond 24 h; treated as exploratory transfer.",
    role == "external" & time_group == "acute_main" ~
      "External acute validation contrast not used as a strict anchor.",
    role == "external" & time_group == "extended_transfer" ~
      "External late-time transfer contrast; exploratory for publication claims.",
    TRUE ~ "Contrast retained for audit with incomplete timing metadata."
  )
}

make_limitation <- function(role, time_group, organism) {
  base <- case_when(
    role == "anchor" ~ "Not independent of the anchor evidence except in leave-one-anchor-out analysis.",
    time_group == "extended_transfer" ~ "Late timepoint; exploratory transfer, not strict acute innate validation.",
    time_group == "unknown_time" ~ "Timing is missing or unclear.",
    role == "calibration" ~ "Calibration dataset; not counted as a strict anchor.",
    TRUE ~ "External dataset; interpret with dataset-specific sample size and design constraints."
  )
  if (!is.na(organism) && nzchar(organism) && organism != "NA") {
    base
  } else {
    paste(base, "Organism metadata unavailable in inspected design files.")
  }
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  eval_tbl <- load_step09_eval(required = TRUE)
  if (nrow(eval_tbl) == 0) stop("Step 09 eval table has no rows.")

  audit_tbl <- map_dfr(seq_len(nrow(eval_tbl)), function(i) {
    row <- eval_tbl[i, ]
    split <- tryCatch(read_split_table(row$split_path), error = function(e) {
      log_msg("WARN", "Could not read split for dataset audit: ", row$split_path,
              " | ", conditionMessage(e))
      tibble()
    })

    organism <- get_any_value(split, c("organism", "species"))
    delivery_platform <- get_value_from_split(
      split, "DELIVERY",
      c("delivery_platform", "delivery_type", "platform", "vector", "formulation",
        "group_raw", "group")
    )
    control_definition <- get_value_from_split(
      split, "CONTROL",
      c("control_definition", "control_label", "group_raw", "group")
    )

    tissue <- if (!is.na(row$tissue) && nzchar(as.character(row$tissue))) {
      as.character(row$tissue)
    } else {
      get_any_value(split, c("tissue", "tissue_raw"))
    }
    time_h <- suppressWarnings(as.numeric(row$time_h))
    role <- assign_audit_role(row$gse_id, time_h)
    time_group <- assign_time_group(time_h)

    tibble(
      gse_id = row$gse_id,
      split_id = row$split_id,
      dataset_type = row$dataset_type,
      organism = organism,
      tissue = tissue,
      time_h = time_h,
      delivery_platform_type = delivery_platform,
      control_definition = control_definition,
      n_controls = suppressWarnings(as.integer(row$n_controls)),
      n_delivery = suppressWarnings(as.integer(row$n_delivery)),
      acute_late_group = time_group,
      anchor_calibration_external = role,
      inclusion_rationale = make_inclusion_rationale(role, time_group),
      limitation = make_limitation(role, time_group, organism),
      split_path = row$split_path
    )
  }) %>%
    arrange(anchor_calibration_external, gse_id, time_h, split_id)

  out_tsv <- file.path(AUDIT_RESULTS_DIR, "dataset_audit_table.tsv")
  out_csv <- file.path(AUDIT_RESULTS_DIR, "dataset_audit_table.csv")
  write_tsv_safe(audit_tbl, out_tsv)
  write_csv_safe(audit_tbl, out_csv)
  log_msg("INFO", "Wrote dataset audit table rows=", nrow(audit_tbl),
          " to ", out_tsv)
  invisible(audit_tbl)
}

main()
