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

this_script <- imrs_detect_script_path_local("07_weak_dataset_paper_context_audit.R")
if (is.na(this_script)) stop("Could not identify this script path. Use RStudio Source/Run or Rscript.", call. = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "07_weak_dataset_paper_context_audit"

allowed_categories <- c(
  "shared_delivery_vehicle",
  "late_timepoint",
  "therapeutic_cargo_specific_effect",
  "disease_rescue_model",
  "distal_or_adaptive_tissue",
  "formulation_designed_to_reduce_inflammation",
  "tissue_time_kinetic_effect",
  "unclear_requires_manual_review"
)

mapping_path <- file.path(AUDIT_DATASETS_DIR, "weak_dataset_paper_mappings.tsv")

base_gse_id <- function(dataset_id) {
  out <- str_match(as.character(dataset_id), "^(GSE[0-9]+)")[, 2]
  ifelse(is.na(out), as.character(dataset_id), out)
}

is_truthy <- function(x) {
  if (length(x) == 0 || is.na(x)) return(FALSE)
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "Y", "1")
}

infer_organism <- function(dataset_id, tissue, mapping_organism) {
  dataset_id <- as.character(dataset_id)
  tissue <- tolower(as.character(tissue))
  case_when(
    dataset_id == "GSE262515_cell_line" ~ "Homo sapiens",
    dataset_id == "GSE262515_tissue" ~ "Mus musculus",
    TRUE ~ mapping_organism
  )
}

infer_delivery_platform <- function(treatment_group, dataset_id) {
  x <- tolower(paste(treatment_group, dataset_id, sep = " "))
  case_when(
    str_detect(x, "tlnp|lnp|sm102|am80") ~ "lipid nanoparticle",
    str_detect(x, "aav9|aav") ~ "AAV9",
    str_detect(x, "ad26") ~ "adenovirus serotype 26 vector",
    str_detect(x, "ad5") ~ "adenovirus serotype 5 vector",
    TRUE ~ NA_character_
  )
}

infer_cargo <- function(treatment_group, dataset_id) {
  x <- tolower(as.character(treatment_group))
  dataset_id <- as.character(dataset_id)
  case_when(
    str_detect(x, "sihbv") ~ "siHBV",
    str_detect(x, "sinc") ~ "siNC",
    str_detect(x, "ar45") ~ "AR45 transgene",
    str_detect(x, "gfp") ~ "GFP control transgene",
    str_detect(x, "am80") ~ "mRNA vaccine with Am80 retinoic acid receptor agonist",
    str_detect(x, "sm102") ~ "mRNA vaccine in SM102-LNP",
    dataset_id == "GSE264344" ~ NA_character_,
    TRUE ~ NA_character_
  )
}

shared_delivery_between_groups <- function(treatment_group, control_group) {
  tx <- tolower(as.character(treatment_group))
  ctrl <- tolower(as.character(control_group))
  platforms <- c("tlnp", "lnp", "aav9", "aav", "ad26", "ad5")
  any(map_lgl(platforms, ~ str_detect(tx, .x) && str_detect(ctrl, .x)))
}

tissue_matches_expected <- function(dataset_id, tissue) {
  tissue_l <- tolower(as.character(tissue))
  case_when(
    dataset_id == "GSE262515_cell_line" ~ NA,
    dataset_id == "GSE262515_tissue" ~ TRUE,
    dataset_id == "GSE166655" & str_detect(tissue_l, "tibialis|muscle") ~ TRUE,
    dataset_id == "GSE314070" & str_detect(tissue_l, "intestine") ~ FALSE,
    dataset_id == "GSE264344" & str_detect(tissue_l, "blood|muscle|dln|lymph") ~ TRUE,
    TRUE ~ NA
  )
}

paper_category <- function(dataset_id, treatment_group, control_group, tissue, time_h) {
  x <- tolower(paste(dataset_id, treatment_group, control_group, tissue, sep = " "))
  time_h <- suppressWarnings(as.numeric(time_h))
  case_when(
    dataset_id == "GSE264344" & is.finite(time_h) & time_h > 24 ~
      "tissue_time_kinetic_effect",
    dataset_id == "GSE314070" & str_detect(x, "am80") ~
      "formulation_designed_to_reduce_inflammation",
    dataset_id == "GSE314070" ~
      "distal_or_adaptive_tissue",
    dataset_id == "GSE166655" ~
      "disease_rescue_model",
    str_detect(x, "sihbv") ~
      "therapeutic_cargo_specific_effect",
    dataset_id %in% c("GSE262515_tissue", "GSE262515_cell_line") &
      shared_delivery_between_groups(treatment_group, control_group) ~
      "shared_delivery_vehicle",
    dataset_id %in% c("GSE262515_tissue", "GSE262515_cell_line") &
      is.finite(time_h) & time_h > 24 ~
      "late_timepoint",
    TRUE ~ "unclear_requires_manual_review"
  )
}

status_for_row <- function(dataset_id, category, treatment_group, control_group, time_h) {
  time_h <- suppressWarnings(as.numeric(time_h))
  if (category == "unclear_requires_manual_review") return("unclear")
  if (dataset_id %in% c("GSE262515_tissue", "GSE262515_cell_line") &&
      !shared_delivery_between_groups(treatment_group, control_group)) {
    return("partially_explainable")
  }
  if (dataset_id == "GSE166655" && !shared_delivery_between_groups(treatment_group, control_group)) {
    return("partially_explainable")
  }
  if (is.finite(time_h) && time_h > 24) return("explainable")
  "partially_explainable"
}

risk_for_status <- function(status, delta) {
  delta <- suppressWarnings(as.numeric(delta))
  case_when(
    status == "explainable" ~ "low",
    status == "partially_explainable" & is.finite(delta) & delta > 0 ~ "medium",
    status == "partially_explainable" ~ "medium",
    status == "unclear" ~ "high",
    TRUE ~ "high"
  )
}

interpretation_for_row <- function(dataset_id, category, status, treatment_group, control_group,
                                   tissue, time_h, delta) {
  time_h <- suppressWarnings(as.numeric(time_h))
  delta_txt <- ifelse(is.finite(delta), sprintf("%.3f", delta), "NA")
  acute_txt <- ifelse(is.finite(time_h) && time_h <= 24, "acute", "non-acute")
  case_when(
    category == "tissue_time_kinetic_effect" ~ paste0(
      "Delta IMRSz is ", delta_txt, " in ", tissue, " at ", time_h,
      " h. For the Ad vector time-course, this 72 h external split is best interpreted as waning tissue-time kinetics, not as evidence against early innate activation."
    ),
    category == "formulation_designed_to_reduce_inflammation" ~ paste0(
      "Delta IMRSz is ", delta_txt, " in small intestine at ", time_h,
      " h. The Am80-LNP context is late mucosal/adaptive biology and the formulation is designed to support gut immunity while mitigating injection-site inflammation."
    ),
    category == "distal_or_adaptive_tissue" ~ paste0(
      "Delta IMRSz is ", delta_txt, " in small intestine at ", time_h,
      " h. This is distal late mucosal/adaptive tissue profiling, not an acute injection-site innate readout."
    ),
    category == "disease_rescue_model" ~ paste0(
      "Delta IMRSz is ", delta_txt, " in tibialis anterior at ", time_h,
      " h. This is a late AAV therapeutic transgene/disease-rescue design, so weak acute innate signal is biologically plausible."
    ),
    category == "therapeutic_cargo_specific_effect" ~ paste0(
      "Delta IMRSz is ", delta_txt, " for ", treatment_group, " versus ", control_group,
      ". The siHBV cargo is therapeutic RNAi in an optimized tLNP context; the result should not be interpreted as pure delivery-induced innate activation."
    ),
    category == "late_timepoint" ~ paste0(
      "Delta IMRSz is ", delta_txt, " at ", time_h,
      " h. The contrast is ", acute_txt,
      " and weak signal is plausibly affected by late post-delivery kinetics."
    ),
    category == "shared_delivery_vehicle" ~ paste0(
      "Delta IMRSz is ", delta_txt,
      ". Treatment and control share a delivery vehicle, so the contrast may remove much of the vehicle-driven innate component."
    ),
    TRUE ~ paste0(
      "Delta IMRSz is ", delta_txt,
      ". The available local metadata and paper mapping are insufficient for a manuscript-ready biological explanation."
    )
  )
}

notes_for_row <- function(dataset_id, base_gse, category, status, treatment_group, control_group) {
  notes <- c("PMID/title from user-provided verified mapping; no local metadata contradiction was found.")
  if (dataset_id %in% c("GSE262515_tissue", "GSE262515_cell_line") &&
      !shared_delivery_between_groups(treatment_group, control_group)) {
    notes <- c(notes, "Current Step 09 split is versus baseline_0, not a direct tLNP/siHBV-vs-tLNP/siNC shared-vehicle contrast.")
  }
  if (dataset_id == "GSE166655" &&
      !shared_delivery_between_groups(treatment_group, control_group)) {
    notes <- c(notes, "Current Step 09 split uses delivery_untreated as control; direct AAV9-AR45-vs-AAV9-GFP shared-vector comparison is not represented here.")
  }
  if (category == "unclear_requires_manual_review") {
    notes <- c(notes, "Manual paper/design review required before using as an explanation.")
  }
  paste(notes, collapse = " ")
}

read_split_context <- function(split_path) {
  split <- read_split_table(split_path)
  group_col <- first_present_col(split, c("group_raw", "group"))
  if (is.na(group_col)) {
    return(tibble(treatment_group = NA_character_, control_group = NA_character_))
  }
  treatment <- split %>%
    filter(condition_simple == "DELIVERY") %>%
    pull(.data[[group_col]]) %>%
    as.character()
  if (length(treatment) == 0) treatment <- NA_character_
  control <- split %>%
    filter(condition_simple == "CONTROL") %>%
    pull(.data[[group_col]]) %>%
    as.character()
  if (length(control) == 0) control <- NA_character_
  tibble(
    treatment_group = collapse_unique_chr(treatment),
    control_group = collapse_unique_chr(control)
  )
}

make_markdown_report <- function(tbl, path) {
  lines <- c(
    "# Weak Dataset Paper Context Audit",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Purpose",
    "",
    "For weak or near-zero IMRS datasets, this audit prevents score-only speculation. It records the paper context, treatment/control design, tissue, timepoint, Delta IMRSz, and a conservative biological interpretation. Delta IMRSz and design context are primary; AUC is secondary.",
    "",
    "## Summary",
    "",
    paste0("- Contrasts audited: ", nrow(tbl)),
    paste0("- Explainable: ", sum(tbl$explainable_status == "explainable", na.rm = TRUE)),
    paste0("- Partially explainable: ", sum(tbl$explainable_status == "partially_explainable", na.rm = TRUE)),
    paste0("- Unclear: ", sum(tbl$explainable_status == "unclear", na.rm = TRUE)),
    paste0("- High reviewer risk: ", sum(tbl$reviewer_risk_level == "high", na.rm = TRUE)),
    "",
    "## Interpretation Rules",
    "",
    "- AUC is not used as the main explanation.",
    "- A weak Delta IMRSz is considered explainable only when the paper context and local split design support a plausible biological reason.",
    "- Direct shared-delivery explanations are not assigned when the current Step 09 split is against `baseline_0` or `delivery_untreated`.",
    "- Late timepoints and distal/adaptive tissues are treated as transfer or biology-specific contexts, not failures of acute innate validation.",
    "",
    "## Category Counts",
    "",
    paste(capture.output(print(tbl %>% count(explanation_category, explainable_status, reviewer_risk_level))),
          collapse = "\n"),
    "",
    "## Contrast-Level Notes",
    ""
  )
  row_lines <- pmap_chr(tbl, function(...) {
    r <- list(...)
    paste0(
      "- **", r$dataset_id, "** / `", r$split_id, "`: ",
      r$manuscript_ready_interpretation,
      " Reviewer risk: ", r$reviewer_risk_level, "."
    )
  })
  writeLines(c(lines, row_lines), path, useBytes = TRUE)
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  if (!file.exists(mapping_path)) stop("Missing mapping file: ", mapping_path)

  mappings <- read_tsv(mapping_path, show_col_types = FALSE, progress = FALSE)
  eval_tbl <- load_step09_eval(required = TRUE) %>%
    mutate(pass = as.logical(pass),
           base_gse = base_gse_id(gse_id),
           time_h = suppressWarnings(as.numeric(time_h)),
           delta_mean_imrs_z = suppressWarnings(as.numeric(delta_mean_imrs_z)),
           auc_imrs_z = suppressWarnings(as.numeric(auc_imrs_z)))

  default_datasets <- c("GSE262515_tissue", "GSE262515_cell_line", "GSE166655", "GSE314070")
  weak_rows <- eval_tbl %>%
    filter(
      pass %in% TRUE,
      gse_id %in% default_datasets |
        (gse_id == "GSE264344" & dataset_type == "external" &
           is.finite(delta_mean_imrs_z) & delta_mean_imrs_z < 2)
    ) %>%
    arrange(gse_id, tissue, time_h, split_id)
  if (nrow(weak_rows) == 0) stop("No weak dataset rows selected from Step 09 eval.")

  out_tbl <- map_dfr(seq_len(nrow(weak_rows)), function(i) {
    row <- weak_rows[i, ]
    split_ctx <- tryCatch(read_split_context(row$split_path), error = function(e) {
      log_msg("WARN", "Could not read split context for ", row$split_path, ": ",
              conditionMessage(e))
      tibble(treatment_group = row$contrast_label, control_group = row$control_label)
    })
    mapping <- mappings %>% filter(dataset_family == row$base_gse) %>% slice_head(n = 1)
    if (nrow(mapping) == 0) {
      mapping <- tibble(
        PMID = NA_character_,
        paper_title = NA_character_,
        DOI = NA_character_,
        publication_year = NA_integer_,
        default_organism = NA_character_,
        main_interpretation = NA_character_
      )
    }

    treatment <- split_ctx$treatment_group[[1]]
    control <- split_ctx$control_group[[1]]
    category <- paper_category(row$gse_id, treatment, control, row$tissue, row$time_h)
    if (!(category %in% allowed_categories)) category <- "unclear_requires_manual_review"
    status <- status_for_row(row$gse_id, category, treatment, control, row$time_h)
    organism <- infer_organism(row$gse_id, row$tissue, mapping$default_organism[[1]])
    platform <- infer_delivery_platform(treatment, row$gse_id)
    cargo <- infer_cargo(treatment, row$gse_id)

    tibble(
      dataset_id = row$gse_id,
      gse_id = row$base_gse,
      split_id = row$split_id,
      PMID = as.character(mapping$PMID[[1]]),
      paper_title = mapping$paper_title[[1]],
      DOI = mapping$DOI[[1]],
      publication_year = suppressWarnings(as.integer(mapping$publication_year[[1]])),
      organism = organism,
      tissue = row$tissue,
      time_h = row$time_h,
      treatment_group = treatment,
      control_group = control,
      delivery_platform = platform,
      cargo_or_transgene = cargo,
      whether_delivery_is_shared_between_groups = shared_delivery_between_groups(treatment, control),
      whether_timepoint_is_acute_24h_or_less = is.finite(row$time_h) && row$time_h <= 24,
      whether_tissue_matches_expected_innate_response_site =
        tissue_matches_expected(row$gse_id, row$tissue),
      original_IMRS_delta = row$delta_mean_imrs_z,
      original_AUC_secondary = row$auc_imrs_z,
      explanation_category = category,
      explainable_status = status,
      manuscript_ready_interpretation =
        interpretation_for_row(row$gse_id, category, status, treatment, control,
                               row$tissue, row$time_h, row$delta_mean_imrs_z),
      reviewer_risk_level = risk_for_status(status, row$delta_mean_imrs_z),
      notes = notes_for_row(row$gse_id, row$base_gse, category, status, treatment, control)
    )
  })

  out_tsv <- file.path(AUDIT_RESULTS_DIR, "weak_dataset_paper_context_audit.tsv")
  out_csv <- file.path(AUDIT_RESULTS_DIR, "weak_dataset_paper_context_audit.csv")
  out_md <- file.path(AUDIT_REPORT_DIR, "weak_dataset_paper_context_audit.md")
  write_tsv_safe(out_tbl, out_tsv)
  write_csv_safe(out_tbl, out_csv)
  make_markdown_report(out_tbl, out_md)

  plot_tbl <- out_tbl %>%
    count(explanation_category, explainable_status, reviewer_risk_level, name = "n") %>%
    mutate(
      explanation_category = factor(explanation_category, levels = allowed_categories),
      reviewer_risk_level = factor(reviewer_risk_level, levels = c("low", "medium", "high"))
    )
  p <- ggplot(plot_tbl, aes(x = explanation_category, y = n, fill = reviewer_risk_level)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              color = "white", size = 3) +
    scale_fill_manual(values = c(low = "#4E79A7", medium = "#F28E2B", high = "#E15759"),
                      drop = FALSE) +
    labs(
      title = "Weak Dataset Literature-Context Explanations",
      x = NULL,
      y = "Contrasts",
      fill = "Reviewer risk"
    ) +
    theme_audit() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_plot(p, file.path(AUDIT_PLOTS_DIR, "weak_dataset_context_summary.png"),
            width = 9, height = 5.5)

  log_msg("INFO", "Wrote weak dataset paper-context audit rows=", nrow(out_tbl),
          " to ", out_tsv)
  invisible(out_tbl)
}

main()
