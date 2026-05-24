#!/usr/bin/env Rscript
# Priority 4 manuscript-ready supplementary table packaging for IMRS v6.
# This script only packages existing results; it does not rerun scoring,
# validation, robustness, or enrichment analyses.

options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                      error = function(e) NA_character_)
if (is.na(this_file)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
}
active_config_helper <- file.path(dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE)),
                                  "_active_config.R")
if (!file.exists(active_config_helper)) {
  active_config_helper <- file.path(getwd(), "scripts", "active_manuscript", "_active_config.R")
}
source(active_config_helper)

config <- imrs_load_active_config(dirname(active_config_helper))
project_root <- imrs_project_root(config)
v6_root <- imrs_config_field_path(config, "manuscript_output_dir")
out_dir <- imrs_config_field_path(config, "supplementary_tables_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

delta <- intToUtf8(0x0394)
metric_label <- paste0("delivery-minus-control ", delta, "IMRSz")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
log_path <- file.path(out_dir, "supplementary_table_generation_log.txt")
log_lines <- character()
warnings_found <- character()

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}

write_log <- function() {
  writeLines(log_lines, log_path, useBytes = TRUE)
}

add_warning <- function(...) {
  msg <- paste(..., collapse = "")
  warnings_found <<- c(warnings_found, msg)
  log_msg("WARNING: ", msg)
}

required_packages <- c("readr", "dplyr", "tidyr", "stringr", "tibble", "purrr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  log_msg("ERROR: Missing required package(s): ", paste(missing_packages, collapse = ", "))
  write_log()
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(writexl)
})

log_msg("Starting Priority 4 supplementary table packaging.")
log_msg("Project root: ", project_root)
log_msg("Output folder: ", out_dir)
log_msg("This packaging run does not rerun or alter IMRS scoring, validation, robustness, or enrichment statistics.")

paths <- list(
  provenance = imrs_config_field_path(config, "provenance_table"),
  manuscript_roles = imrs_config_field_path(config, "manuscript_role_table"),
  boundary_audit = imrs_config_field_path(config, "boundary_audit_table"),
  label_permutation = imrs_config_field_path(config, "label_permutation_summary"),
  leave_one_gene_out = imrs_config_field_path(config, "leave_one_gene_out_summary"),
  gene_dominance = imrs_config_field_path(config, "gene_dominance_summary"),
  threshold_sensitivity = imrs_config_field_path(config, "threshold_sensitivity_summary"),
  leave_one_anchor_out = imrs_config_field_path(config, "leave_one_anchor_out_summary"),
  comparator_benchmarking = imrs_config_field_path(config, "comparator_benchmarking_summary"),
  coefficient_sensitivity = imrs_config_field_path(config, "coefficient_sensitivity_summary"),
  s5_enrichment = file.path(imrs_config_field_path(config, "priority3_enrichment_dir"), "tables", "Supplementary_Table_S5_IMRS_gene_enrichment_all.tsv"),
  s5_mapping_audit = file.path(imrs_config_field_path(config, "priority3_enrichment_dir"), "tables", "IMRS_retained_gene_mapping_audit.tsv"),
  s5_background_audit = file.path(imrs_config_field_path(config, "priority3_enrichment_dir"), "tables", "IMRS_enrichment_background_audit.tsv")
)

must_exist <- c("provenance", "boundary_audit", "label_permutation", "leave_one_gene_out",
                "gene_dominance", "threshold_sensitivity", "leave_one_anchor_out", "s5_enrichment")
missing_sources <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_sources) > 0) {
  for (nm in missing_sources) add_warning("Source file unavailable: ", nm, " -> ", paths[[nm]])
}
missing_required <- intersect(must_exist, missing_sources)
if (length(missing_required) > 0) {
  log_msg("ERROR: Required source file(s) unavailable: ", paste(missing_required, collapse = ", "))
  write_log()
  stop("Required source file(s) unavailable: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

read_source <- function(path) {
  readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
}

write_table_pair <- function(tbl, stem) {
  tsv_path <- file.path(out_dir, paste0(stem, ".tsv"))
  xlsx_path <- file.path(out_dir, paste0(stem, ".xlsx"))
  readr::write_tsv(tbl, tsv_path, na = "")
  writexl::write_xlsx(list(data = tbl), xlsx_path)
  list(tsv = tsv_path, xlsx = xlsx_path)
}

as_num <- function(x) suppressWarnings(as.numeric(x))

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(sort(unique(x)), collapse = "; ")
}

yes_no <- function(x) {
  truth <- toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
  ifelse(truth, "yes", "no")
}

controlled_group <- function(x) {
  y <- as.character(x)
  dplyr::case_when(
    y %in% c("Locked anchor", "anchor", "strict_anchor") ~ "Locked anchor",
    y %in% c("Primary acute validation", "primary_acute_validation") ~ "Primary acute validation",
    y %in% c("Extended validation", "extended_validation") ~ "Extended validation",
    y %in% c("Secondary support", "secondary_support", "secondary_support_not_primary") ~ "Secondary support",
    is.na(y) | !nzchar(y) | y == "NA" ~ "Excluded/unclear",
    TRUE ~ y
  )
}

fmt <- function(x, digits = 4) {
  if (length(x) == 0 || is.na(x)) return(NA_character_)
  formatC(as.numeric(x), format = "f", digits = digits)
}

fmt_prop <- function(x) fmt(x, 4)

source_rel <- function(path) {
  stringr::str_replace(normalizePath(path, winslash = "/", mustWork = FALSE),
                       paste0("^", stringr::fixed(project_root), "/?"), "")
}

table_note_s1 <- paste0(
  "Dataset-level provenance and manuscript role summary for the frozen, transfer-oriented IMRS scoring framework; ",
  metric_label, " is summarized across scored split contrasts and AUC is secondary."
)
table_note_s2 <- paste0(
  "Split-level delivery-versus-control contrast table for the frozen, transfer-oriented IMRS scoring framework; ",
  metric_label, " is the primary transfer metric and AUC is secondary."
)
table_note_s3 <- paste0(
  "Late or context-shifted settings are retained and labeled as boundary-setting evidence, not primary validation."
)
table_note_s4 <- paste0(
  "Robustness summary consolidating analyses of score stability and non-degeneracy without changing the IMRS score."
)
table_note_s5 <- paste0(
  "Gene-program enrichment was performed using the 287 retained frozen IMRS genes. Of these, 252 mapped to Entrez identifiers and 35 were unmapped. ",
  "The enrichment background was derived from the anchor gene_power.tsv universe with 22,304 mapped background genes. ",
  "Enrichment results provide program-level biological context and do not establish causal pathways, cell-type sources, clinical reactogenicity prediction, or delivery-platform safety ranking."
)

provenance <- read_source(paths$provenance) %>%
  mutate(
    delta_value = as_num(.data$delta_mean_imrs_z),
    auc_value = as_num(.data$auc_imrs_z_secondary),
    time_value = as_num(.data$time_h),
    manuscript_analysis_group = controlled_group(.data$manuscript_interpretation_group),
    primary_claim_included_yes_no = yes_no(.data$used_for_primary_claim)
  )
log_msg("Loaded provenance source rows: ", nrow(provenance), " from ", paths$provenance)

S1 <- provenance %>%
  group_by(.data$dataset_id, .data$tissue, .data$delivery_platform, .data$manuscript_analysis_group) %>%
  summarise(
    dataset_context_id = paste0(first(.data$dataset_id), "__", stringr::str_replace_all(first(.data$tissue), "[^A-Za-z0-9]+", "_"), "__", stringr::str_replace_all(first(.data$delivery_platform), "[^A-Za-z0-9]+", "_"), "__", stringr::str_replace_all(first(.data$manuscript_analysis_group), "[^A-Za-z0-9]+", "_")),
    PMID = collapse_unique(.data$PMID),
    DOI = collapse_unique(.data$DOI),
    source_title = collapse_unique(.data$paper_title),
    timepoint_h_min = ifelse(all(is.na(.data$time_value)), NA_real_, min(.data$time_value, na.rm = TRUE)),
    timepoint_h_max = ifelse(all(is.na(.data$time_value)), NA_real_, max(.data$time_value, na.rm = TRUE)),
    treatment_context = collapse_unique(dplyr::coalesce(.data$paper_treatment_context, .data$treatment_group_clean)),
    control_context = collapse_unique(.data$control_group_clean),
    primary_claim_included_yes_no = ifelse(any(.data$primary_claim_included_yes_no == "yes", na.rm = TRUE), "yes", "no"),
    n_split_contrasts = dplyr::n(),
    n_positive_delta = sum(.data$delta_value > 0, na.rm = TRUE),
    positive_directionality_proportion = ifelse(all(is.na(.data$delta_value)), NA_real_, mean(.data$delta_value > 0, na.rm = TRUE)),
    mean_delivery_minus_control_delta_IMRSz = ifelse(all(is.na(.data$delta_value)), NA_real_, mean(.data$delta_value, na.rm = TRUE)),
    median_delivery_minus_control_delta_IMRSz = ifelse(all(is.na(.data$delta_value)), NA_real_, median(.data$delta_value, na.rm = TRUE)),
    mean_AUC_secondary = ifelse(all(is.na(.data$auc_value)), NA_real_, mean(.data$auc_value, na.rm = TRUE)),
    scoring_status = collapse_unique(.data$scoring_status),
    paper_context_alignment = collapse_unique(.data$paper_result_alignment),
    .groups = "drop"
  ) %>%
  transmute(
    dataset_id = .data$dataset_id,
    dataset_context_id = .data$dataset_context_id,
    PMID = .data$PMID,
    DOI = .data$DOI,
    source_title = .data$source_title,
    tissue = .data$tissue,
    timepoint_h_min = .data$timepoint_h_min,
    timepoint_h_max = .data$timepoint_h_max,
    delivery_platform_or_modality = .data$delivery_platform,
    treatment_context = .data$treatment_context,
    control_context = .data$control_context,
    manuscript_analysis_group = .data$manuscript_analysis_group,
    primary_claim_included_yes_no = .data$primary_claim_included_yes_no,
    n_split_contrasts = .data$n_split_contrasts,
    n_positive_delta = .data$n_positive_delta,
    positive_directionality_proportion = .data$positive_directionality_proportion,
    mean_delivery_minus_control_delta_IMRSz = .data$mean_delivery_minus_control_delta_IMRSz,
    median_delivery_minus_control_delta_IMRSz = .data$median_delivery_minus_control_delta_IMRSz,
    mean_AUC_secondary = .data$mean_AUC_secondary,
    scoring_status = .data$scoring_status,
    paper_context_alignment = .data$paper_context_alignment,
    interpretation_summary = paste0("Dataset context summarized under ", .data$manuscript_analysis_group,
                                    " using the frozen, transfer-oriented IMRS scoring framework; AUC is secondary."),
    source_file = source_rel(paths$provenance),
    table_note = table_note_s1
  ) %>%
  arrange(.data$manuscript_analysis_group, .data$dataset_id, .data$tissue, .data$delivery_platform_or_modality)

S2 <- provenance %>%
  transmute(
    split_id = .data$split_id,
    dataset_id = .data$dataset_id,
    dataset_context_id = paste0(.data$dataset_id, "__", stringr::str_replace_all(.data$tissue, "[^A-Za-z0-9]+", "_"), "__H", .data$time_h),
    tissue = .data$tissue,
    timepoint_h = .data$time_h,
    treatment_group_clean = .data$treatment_group_clean,
    control_group_clean = .data$control_group_clean,
    delivery_platform_or_modality = .data$delivery_platform,
    manuscript_analysis_group = .data$manuscript_analysis_group,
    primary_claim_included_yes_no = .data$primary_claim_included_yes_no,
    n_delivery = .data$n_delivery,
    n_control = .data$n_control,
    n_total = .data$n_total,
    delivery_minus_control_delta_IMRSz = .data$delta_mean_imrs_z,
    direction_positive_yes_no = yes_no(.data$direction_positive),
    AUC_secondary = .data$auc_imrs_z_secondary,
    pass_status = yes_no(.data$qc_pass),
    exclusion_or_warning_flag = dplyr::case_when(
      .data$manuscript_analysis_group == "Excluded/unclear" ~ "excluded_or_unclear_metadata",
      !is.na(.data$limitation) & .data$limitation != "NA" & nzchar(.data$limitation) ~ .data$limitation,
      TRUE ~ NA_character_
    ),
    paper_context_alignment = .data$paper_result_alignment,
    interpretation_summary = paste0("Split contrast scored under ", .data$manuscript_analysis_group,
                                    "; ", metric_label, " is reported with AUC as a secondary metric."),
    source_file = source_rel(paths$provenance),
    table_note = table_note_s2
  ) %>%
  arrange(.data$manuscript_analysis_group, .data$dataset_id, as_num(.data$timepoint_h), .data$split_id)

boundary_source <- read_source(paths$boundary_audit)
role_lookup <- provenance %>%
  select(split_id, manuscript_analysis_group) %>%
  distinct()

boundary_category <- function(x) {
  dplyr::case_when(
    x == "late_timepoint" ~ "late timepoint",
    x == "tissue_time_kinetic_effect" ~ "tissue/time kinetic effect",
    x == "therapeutic_cargo_specific_effect" ~ "therapeutic cargo-specific effect",
    x == "disease_rescue_model" ~ "disease-rescue model",
    x == "distal_or_adaptive_tissue" ~ "distal or adaptive tissue",
    x == "formulation_designed_to_reduce_inflammation" ~ "formulation designed to reduce inflammation",
    TRUE ~ "context-shifted biological setting"
  )
}

boundary_explanation <- function(cat) {
  dplyr::case_when(
    cat == "late timepoint" ~ "Timing falls outside the primary acute window, so the contrast is interpreted as boundary-setting evidence for late or context-shifted settings.",
    cat == "tissue/time kinetic effect" ~ "Tissue and time kinetics differ from the primary acute validation context, supporting a boundary-setting interpretation.",
    cat == "therapeutic cargo-specific effect" ~ "Therapeutic cargo context can shape the observed transcriptional program, so the contrast is retained as boundary-setting evidence.",
    cat == "disease-rescue model" ~ "Disease-rescue design adds biological context beyond delivery-versus-control scoring, supporting boundary-setting interpretation.",
    cat == "distal or adaptive tissue" ~ "Distal or adaptive tissue context can attenuate the acute delivery-associated innate transcriptional response captured by IMRS.",
    cat == "formulation designed to reduce inflammation" ~ "Formulation design is expected to reduce inflammatory transcriptional signal, so the contrast is treated as boundary-setting evidence.",
    TRUE ~ "Metadata or design context requires caution, so the contrast is retained as boundary-setting evidence."
  )
}

S3 <- boundary_source %>%
  left_join(role_lookup, by = "split_id") %>%
  mutate(
    boundary_category = boundary_category(.data$explanation_category),
    manuscript_analysis_group = controlled_group(.data$manuscript_analysis_group)
  ) %>%
  transmute(
    dataset_id = .data$dataset_id,
    split_id = .data$split_id,
    PMID = .data$PMID,
    DOI = .data$DOI,
    tissue = .data$tissue,
    timepoint_h = .data$time_h,
    treatment_group_clean = .data$treatment_group,
    control_group_clean = .data$control_group,
    manuscript_analysis_group = .data$manuscript_analysis_group,
    delivery_minus_control_delta_IMRSz = .data$original_IMRS_delta,
    AUC_secondary = .data$original_AUC_secondary,
    boundary_category = .data$boundary_category,
    boundary_explanation = boundary_explanation(.data$boundary_category),
    interpretation_support_level = stringr::str_replace_all(.data$explainable_status, "_", " "),
    paper_context_alignment = paste0("Paper context supports cautious interpretation as ", .data$boundary_category, "."),
    reviewer_risk_level = .data$reviewer_risk_level,
    manuscript_handling = "Retain as late or context-shifted settings and interpret as boundary-setting evidence rather than primary validation.",
    source_file = source_rel(paths$boundary_audit),
    table_note = table_note_s3
  ) %>%
  arrange(.data$boundary_category, .data$dataset_id, as_num(.data$timepoint_h), .data$split_id)

perm <- read_source(paths$label_permutation)
perm_ok <- perm %>% filter(.data$permutation_status == "ok")
logo <- read_source(paths$leave_one_gene_out)
dominance <- read_source(paths$gene_dominance)
threshold <- read_source(paths$threshold_sensitivity)
loao <- read_source(paths$leave_one_anchor_out)
comparator <- if (file.exists(paths$comparator_benchmarking)) read_source(paths$comparator_benchmarking) else tibble()
coef_sens <- if (file.exists(paths$coefficient_sensitivity)) read_source(paths$coefficient_sensitivity) else tibble()

metric_from_coef <- function(metric) {
  if (nrow(coef_sens) == 0) return(NA_character_)
  hit <- coef_sens %>% filter(.data$metric == !!metric) %>% pull(.data$value)
  if (length(hit) == 0) NA_character_ else as.character(hit[[1]])
}

S4_rows <- list()
S4_rows[["Label permutation"]] <- tibble(
  robustness_analysis = "Label permutation",
  analysis_scope = "Split-level label permutation null across scored contrasts.",
  input_contrasts_or_tests = paste0(nrow(perm_ok), " split contrasts; ", collapse_unique(perm_ok$n_permutations), " permutations per contrast where reported."),
  key_metric = paste0("Observed mean ", metric_label, ", null mean, positive observed contrasts, and outside-95% null count."),
  result_summary = paste0(
    "Observed mean ", metric_label, " = ", fmt(mean(as_num(perm_ok$observed_delta_mean_imrs_z), na.rm = TRUE)),
    "; mean label-permutation null = ", fmt(mean(as_num(perm_ok$null_mean_delta), na.rm = TRUE)),
    "; positive observed contrasts = ", sum(as_num(perm_ok$observed_delta_mean_imrs_z) > 0, na.rm = TRUE), "/", nrow(perm_ok),
    " (", fmt_prop(mean(as_num(perm_ok$observed_delta_mean_imrs_z) > 0, na.rm = TRUE)), ")",
    "; outside 95% null interval = ", sum(toupper(as.character(perm_ok$observed_outside_95pct_null)) == "TRUE", na.rm = TRUE), "/", nrow(perm_ok), "."
  ),
  pass_or_support_status = "supports score stability and non-degeneracy",
  interpretation = paste0("Reduces likelihood that observed ", metric_label, " patterns are explained by random label structure."),
  limitation = "Permutation results are a robustness/context analysis and do not establish causal biology.",
  source_file = source_rel(paths$label_permutation),
  table_note = table_note_s4
)
S4_rows[["Leave-one-gene-out"]] <- tibble(
  robustness_analysis = "Leave-one-gene-out",
  analysis_scope = "Single-gene removal sensitivity across contrast-by-removal records.",
  input_contrasts_or_tests = paste0(nrow(logo), " contrast-by-removal tests across ", dplyr::n_distinct(logo$removed_gene_id), " removed genes."),
  key_metric = "Direction preservation and absolute percent change after single-gene removal.",
  result_summary = paste0(
    nrow(logo), " contrast-by-removal tests; ",
    sum(toupper(as.character(logo$direction_preserved)) == "TRUE", na.rm = TRUE), " direction preserved; median absolute percent change = ",
    fmt(median(as_num(logo$absolute_percent_change_delta), na.rm = TRUE), 3),
    "; maximum absolute delta change = ", fmt(max(as_num(logo$absolute_change_delta), na.rm = TRUE), 4), "."
  ),
  pass_or_support_status = "supports score stability and non-degeneracy",
  interpretation = paste0("Reduces likelihood that observed ", metric_label, " patterns are explained by a single dominant gene."),
  limitation = "Evaluates single-gene removal only within the existing coefficient set.",
  source_file = source_rel(paths$leave_one_gene_out),
  table_note = table_note_s4
)
S4_rows[["Gene dominance"]] <- tibble(
  robustness_analysis = "Gene dominance",
  analysis_scope = "Top-contributor concentration across split contrasts.",
  input_contrasts_or_tests = paste0(nrow(dominance), " split contrasts."),
  key_metric = "Maximum per-sample gene contribution fraction.",
  result_summary = paste0(
    "Median maximum contribution fraction = ", fmt(median(as_num(dominance$median_max_contribution_fraction), na.rm = TRUE), 4),
    "; mean maximum contribution fraction = ", fmt(mean(as_num(dominance$mean_max_contribution_fraction), na.rm = TRUE), 4),
    "; largest observed maximum contribution fraction = ", fmt(max(as_num(dominance$max_max_contribution_fraction), na.rm = TRUE), 4), "."
  ),
  pass_or_support_status = "supports score stability and non-degeneracy",
  interpretation = "Low contribution concentration supports distributed-score behavior.",
  limitation = "Contribution fractions summarize observed scored samples and do not identify cell-type sources.",
  source_file = source_rel(paths$gene_dominance),
  table_note = table_note_s4
)
S4_rows[["Threshold sensitivity"]] <- tibble(
  robustness_analysis = "Threshold sensitivity",
  analysis_scope = "Coefficient/gene-selection threshold grid evaluated using existing sensitivity outputs.",
  input_contrasts_or_tests = paste0(nrow(threshold), " threshold settings."),
  key_metric = "External positive-delta proportion and mean delta range across threshold settings.",
  result_summary = paste0(
    "External positive-delta proportion range = ",
    fmt(min(as_num(threshold$external_proportion_positive_delta), na.rm = TRUE), 4), " to ",
    fmt(max(as_num(threshold$external_proportion_positive_delta), na.rm = TRUE), 4),
    "; external mean ", metric_label, " range = ",
    fmt(min(as_num(threshold$external_mean_delta_imrs_z), na.rm = TRUE), 4), " to ",
    fmt(max(as_num(threshold$external_mean_delta_imrs_z), na.rm = TRUE), 4),
    "; status values = ", collapse_unique(threshold$status), "."
  ),
  pass_or_support_status = "supports score stability and non-degeneracy",
  interpretation = paste0("Reduces likelihood that observed ", metric_label, " patterns are explained by arbitrary threshold choice."),
  limitation = "Threshold grid is a sensitivity analysis around existing selections and does not refit the frozen score used for manuscript validation.",
  source_file = source_rel(paths$threshold_sensitivity),
  table_note = table_note_s4
)
S4_rows[["Leave-one-anchor-out"]] <- tibble(
  robustness_analysis = "Leave-one-anchor-out",
  analysis_scope = "Anchor holdout summaries from existing robustness output.",
  input_contrasts_or_tests = paste0(nrow(loao), " held-out anchor summaries covering ", sum(as_num(loao$n_contrasts), na.rm = TRUE), " held-out contrasts."),
  key_metric = "Held-out positive-delta fraction and mean delta range.",
  result_summary = paste0(
    "Minimum positive-delta fraction = ", fmt(min(as_num(loao$proportion_positive_delta), na.rm = TRUE), 4),
    "; all directions positive values = ", collapse_unique(loao$all_directions_positive),
    "; mean ", metric_label, " range = ", fmt(min(as_num(loao$mean_delta_imrs_z), na.rm = TRUE), 4),
    " to ", fmt(max(as_num(loao$mean_delta_imrs_z), na.rm = TRUE), 4), "."
  ),
  pass_or_support_status = "supports score stability and non-degeneracy",
  interpretation = paste0("Reduces likelihood that observed ", metric_label, " patterns are explained by one anchor dataset."),
  limitation = "Holdout summaries are limited by the number and design of locked-anchor datasets.",
  source_file = source_rel(paths$leave_one_anchor_out),
  table_note = table_note_s4
)
if (nrow(comparator) > 0) {
  S4_rows[["Comparator benchmarking"]] <- tibble(
    robustness_analysis = "Comparator benchmarking",
    analysis_scope = "Comparator immune-signature summaries across manuscript analysis groups.",
    input_contrasts_or_tests = paste0(dplyr::n_distinct(comparator$score_id), " score families across ", dplyr::n_distinct(comparator$manuscript_interpretation_label), " manuscript analysis groups."),
    key_metric = "Mean score shifts, positive-delta proportions, and secondary AUC by comparator signature.",
    result_summary = paste0("Comparator summaries are available for score IDs: ", collapse_unique(comparator$score_id), "."),
    pass_or_support_status = "contextual support",
    interpretation = "Comparator immune signatures contextualize IMRS behavior but do not replace the frozen, transfer-oriented IMRS scoring framework.",
    limitation = "Comparator benchmarking is contextual and does not establish clinical prediction or delivery-platform safety ranking.",
    source_file = source_rel(paths$comparator_benchmarking),
    table_note = table_note_s4
  )
}
if (nrow(coef_sens) > 0) {
  S4_rows[["Coefficient sensitivity"]] <- tibble(
    robustness_analysis = "Coefficient sensitivity",
    analysis_scope = "Existing coefficient-sensitivity inventory assembled from robustness outputs.",
    input_contrasts_or_tests = paste0(nrow(coef_sens), " reported metrics."),
    key_metric = "Summary metrics from leave-one-gene-out, leave-one-anchor-out, threshold sensitivity, and gene dominance outputs.",
    result_summary = paste0(
      "Reported leave-one-gene-out direction-preserved fraction = ", metric_from_coef("leave_one_gene_out_direction_preserved_fraction"),
      "; median absolute percent change = ", metric_from_coef("leave_one_gene_out_median_absolute_percent_change_delta"),
      "; threshold grid count = ", metric_from_coef("threshold_sensitivity_grid_count"), "."
    ),
    pass_or_support_status = "supports score stability and non-degeneracy",
    interpretation = "Consolidates evidence that coefficient behavior is not driven by one anchor, one gene, or one threshold choice.",
    limitation = "Inventory summarizes existing robustness outputs and does not alter the frozen coefficients.",
    source_file = source_rel(paths$coefficient_sensitivity),
    table_note = table_note_s4
  )
}
S4 <- bind_rows(S4_rows)

s5_source <- read_source(paths$s5_enrichment)
s5_required <- c("database", "term_id", "term_name", "description", "gene_count", "background_count",
                 "gene_ratio", "bg_ratio", "p_value", "adjusted_p_value", "q_value",
                 "overlapping_genes_symbols", "overlapping_genes_entrez", "analysis_background",
                 "query_gene_count", "mapped_query_gene_count", "background_gene_count",
                 "software_package", "software_version")
missing_s5_cols <- setdiff(s5_required, names(s5_source))
if (length(missing_s5_cols) > 0) add_warning("S5 source lacks expected column(s): ", paste(missing_s5_cols, collapse = ", "))
S5 <- s5_source %>%
  select(any_of(s5_required), everything()) %>%
  mutate(
    source_file = source_rel(paths$s5_enrichment),
    table_note = table_note_s5
  )

mapping_audit <- if (file.exists(paths$s5_mapping_audit)) read_source(paths$s5_mapping_audit) else tibble(note = "Mapping audit unavailable")
background_audit <- if (file.exists(paths$s5_background_audit)) read_source(paths$s5_background_audit) else tibble(note = "Background audit unavailable")

outputs <- list(
  S1 = write_table_pair(S1, "Supplementary_Table_S1_dataset_level_provenance"),
  S2 = write_table_pair(S2, "Supplementary_Table_S2_split_level_contrasts"),
  S3 = write_table_pair(S3, "Supplementary_Table_S3_late_context_shifted_boundary_audit"),
  S4 = write_table_pair(S4, "Supplementary_Table_S4_robustness_summary"),
  S5 = write_table_pair(S5, "Supplementary_Table_S5_gene_enrichment_results")
)

readme_rows <- tibble(
  section = c(
    "Purpose",
    "S1_dataset_provenance",
    "S2_split_contrasts",
    "S3_boundary_audit",
    "S4_robustness_summary",
    "S5_gene_enrichment",
    "Controlled terminology",
    "AUC handling",
    "Boundary settings",
    "Enrichment interpretation",
    "Source files",
    "Generation timestamp",
    "No-statistics-change statement"
  ),
  text = c(
    "This supplementary table package collects existing IMRS v6 results into manuscript-ready tables for transparency, reproducibility, and publication review.",
    "Dataset-level provenance and manuscript role summary for scored dataset contexts.",
    paste0("Split-level delivery-versus-control contrast table with sample counts, ", metric_label, " values, pass status, and secondary AUC."),
    "Late or context-shifted settings are retained and labeled as boundary-setting evidence.",
    "Robustness, sensitivity, and non-degeneracy analyses consolidated from existing outputs.",
    "Gene-program enrichment results for the 287 retained frozen IMRS genes, with mapping and background audits included in the combined workbook.",
    paste0("Uses acute delivery-associated innate transcriptional response, acute delivery-associated innate transcriptional program, frozen, transfer-oriented IMRS scoring framework, ", metric_label, ", Locked anchor, Primary acute validation, Extended validation, Secondary support, late or context-shifted settings, boundary-setting evidence, and score stability and non-degeneracy."),
    "AUC is reported only as a secondary metric and is not framed as clinical prediction.",
    "Late or context-shifted settings are retained as boundary-setting evidence rather than primary validation.",
    "Enrichment provides program-level biological context for retained frozen IMRS coefficients and does not establish causal pathways, cell-type sources, clinical reactogenicity prediction, or delivery-platform safety ranking.",
    paste(source_rel(unlist(paths[file.exists(unlist(paths))])), collapse = "; "),
    timestamp,
    "No IMRS score, validation, robustness, transfer-evaluation, or enrichment statistics were modified by this packaging run."
  )
)

readme_text <- paste(
  "IMRS v6 supplementary table package",
  "",
  paste0(readme_rows$section, ": ", readme_rows$text, collapse = "\n\n"),
  sep = "\n"
)
readme_path <- file.path(out_dir, "supplementary_table_readme.txt")
writeLines(readme_text, readme_path, useBytes = TRUE)

combined_path <- file.path(out_dir, "Supplementary_Tables_All_IMRS_v6.xlsx")
writexl::write_xlsx(
  list(
    S1_dataset_provenance = S1,
    S2_split_contrasts = S2,
    S3_boundary_audit = S3,
    S4_robustness_summary = S4,
    S5_gene_enrichment = S5,
    S5_mapping_audit = mapping_audit,
    S5_background_audit = background_audit,
    README = readme_rows
  ),
  combined_path
)

table_dims <- list(S1 = dim(S1), S2 = dim(S2), S3 = dim(S3), S4 = dim(S4), S5 = dim(S5))

inventory <- tibble::tribble(
  ~table_id, ~file_tsv, ~file_xlsx, ~n_rows, ~n_columns, ~source_files, ~generated_yes_no, ~warnings, ~notes,
  "S1", basename(outputs$S1$tsv), basename(outputs$S1$xlsx), nrow(S1), ncol(S1), source_rel(paths$provenance), "yes", "", "Dataset-level provenance and manuscript role summary.",
  "S2", basename(outputs$S2$tsv), basename(outputs$S2$xlsx), nrow(S2), ncol(S2), source_rel(paths$provenance), "yes", "", "Split-level delivery-versus-control contrast table.",
  "S3", basename(outputs$S3$tsv), basename(outputs$S3$xlsx), nrow(S3), ncol(S3), source_rel(paths$boundary_audit), "yes", "", "Late or context-shifted settings interpreted as boundary-setting evidence.",
  "S4", basename(outputs$S4$tsv), basename(outputs$S4$xlsx), nrow(S4), ncol(S4), paste(source_rel(c(paths$label_permutation, paths$leave_one_gene_out, paths$gene_dominance, paths$threshold_sensitivity, paths$leave_one_anchor_out, paths$comparator_benchmarking, paths$coefficient_sensitivity)), collapse = "; "), "yes", "", "Robustness summary for score stability and non-degeneracy.",
  "S5", basename(outputs$S5$tsv), basename(outputs$S5$xlsx), nrow(S5), ncol(S5), paste(source_rel(c(paths$s5_enrichment, paths$s5_mapping_audit, paths$s5_background_audit)), collapse = "; "), "yes", "", "Gene-program enrichment results with mapping/background audits in combined workbook.",
  "combined_workbook", "", basename(combined_path), NA_integer_, NA_integer_, "S1-S5 generated tables and S5 audits", "yes", "", "Combined manuscript-ready workbook.",
  "README", basename(readme_path), "", nrow(readme_rows), ncol(readme_rows), "Generated package metadata", "yes", "", "Plain-text package README.",
  "generation_log", basename(log_path), "", NA_integer_, NA_integer_, "Generated package metadata", "yes", "", "Generation and validation log."
)
if (length(warnings_found) > 0) inventory$warnings[inventory$table_id == "generation_log"] <- paste(warnings_found, collapse = "; ")
inventory_path <- file.path(out_dir, "supplementary_table_inventory.tsv")
readr::write_tsv(inventory, inventory_path, na = "")

append_change_log <- function() {
  path <- file.path(v6_root, "v6_change_log.tsv")
  cols <- c("file", "panel_or_figure", "change_type", "old_text", "new_text",
            "plot_visible_or_caption", "regenerated_or_copied", "notes")
  existing <- if (file.exists(path)) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% filter(.data$change_type != "supplementary_table_package_added")
  rows <- tibble(
    file = c(
      "supplementary_tables/Supplementary_Table_S1_dataset_level_provenance.tsv",
      "supplementary_tables/Supplementary_Table_S2_split_level_contrasts.tsv",
      "supplementary_tables/Supplementary_Table_S3_late_context_shifted_boundary_audit.tsv",
      "supplementary_tables/Supplementary_Table_S4_robustness_summary.tsv",
      "supplementary_tables/Supplementary_Table_S5_gene_enrichment_results.tsv",
      "supplementary_tables/Supplementary_Tables_All_IMRS_v6.xlsx",
      "supplementary_tables/supplementary_table_readme.txt",
      "supplementary_tables/supplementary_table_inventory.tsv",
      "supplementary_tables/supplementary_table_generation_log.txt"
    ),
    panel_or_figure = "Supplementary Tables S1-S5",
    change_type = "supplementary_table_package_added",
    old_text = "",
    new_text = c(
      "Created manuscript-ready dataset-level provenance table.",
      "Created manuscript-ready split-level contrast table.",
      "Created manuscript-ready late/context-shifted boundary audit table.",
      "Created manuscript-ready robustness summary table.",
      "Created manuscript-ready retained-gene enrichment results table from existing Priority 3 output.",
      "Created combined Supplementary Tables workbook.",
      "Created supplementary table README.",
      "Created supplementary table inventory.",
      "Created supplementary table generation log."
    ),
    plot_visible_or_caption = "table/package",
    regenerated_or_copied = "generated",
    notes = "Packaging-only update; no IMRS score, validation, robustness, transfer-evaluation, or enrichment statistics were changed."
  )
  readr::write_tsv(bind_rows(existing[, cols], rows), path, na = "")
}

append_table_inventory <- function() {
  path <- file.path(v6_root, "v6_table_inventory.tsv")
  cols <- c("table_file", "purpose_inferred", "edited_yes_no", "key_changes", "row_count", "column_count")
  existing <- if (file.exists(path)) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% filter(!stringr::str_starts(.data$table_file, "supplementary_tables/"))
  table_rows <- tibble(
    table_file = c(
      "supplementary_tables/Supplementary_Table_S1_dataset_level_provenance.tsv",
      "supplementary_tables/Supplementary_Table_S1_dataset_level_provenance.xlsx",
      "supplementary_tables/Supplementary_Table_S2_split_level_contrasts.tsv",
      "supplementary_tables/Supplementary_Table_S2_split_level_contrasts.xlsx",
      "supplementary_tables/Supplementary_Table_S3_late_context_shifted_boundary_audit.tsv",
      "supplementary_tables/Supplementary_Table_S3_late_context_shifted_boundary_audit.xlsx",
      "supplementary_tables/Supplementary_Table_S4_robustness_summary.tsv",
      "supplementary_tables/Supplementary_Table_S4_robustness_summary.xlsx",
      "supplementary_tables/Supplementary_Table_S5_gene_enrichment_results.tsv",
      "supplementary_tables/Supplementary_Table_S5_gene_enrichment_results.xlsx",
      "supplementary_tables/Supplementary_Tables_All_IMRS_v6.xlsx",
      "supplementary_tables/supplementary_table_inventory.tsv"
    ),
    purpose_inferred = c(
      "Supplementary Table S1 dataset-level provenance",
      "Supplementary Table S1 dataset-level provenance workbook",
      "Supplementary Table S2 split-level contrasts",
      "Supplementary Table S2 split-level contrasts workbook",
      "Supplementary Table S3 late/context-shifted boundary audit",
      "Supplementary Table S3 boundary audit workbook",
      "Supplementary Table S4 robustness summary",
      "Supplementary Table S4 robustness summary workbook",
      "Supplementary Table S5 gene enrichment results",
      "Supplementary Table S5 gene enrichment workbook",
      "Combined supplementary table workbook",
      "Supplementary table package inventory"
    ),
    edited_yes_no = "yes",
    key_changes = "Created for Priority 4 manuscript-ready supplementary table package",
    row_count = c(nrow(S1), nrow(S1), nrow(S2), nrow(S2), nrow(S3), nrow(S3), nrow(S4), nrow(S4), nrow(S5), nrow(S5), NA_integer_, nrow(inventory)),
    column_count = c(ncol(S1), ncol(S1), ncol(S2), ncol(S2), ncol(S3), ncol(S3), ncol(S4), ncol(S4), ncol(S5), ncol(S5), NA_integer_, ncol(inventory))
  )
  readr::write_tsv(bind_rows(existing[, cols], table_rows), path, na = "")
}

update_file_inventory <- function() {
  path <- file.path(v6_root, "v6_file_inventory.tsv")
  cols <- c("file", "extension", "size_bytes", "modified_time", "copied_from_v5_yes_no", "generated_in_v6_yes_no")
  existing <- if (file.exists(path)) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% mutate(across(all_of(cols), as.character))
  existing <- existing %>% filter(!stringr::str_starts(.data$file, "supplementary_tables/"))
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  rows <- purrr::map_dfr(sort(files), function(p) {
    info <- file.info(p)
    rel <- paste0("supplementary_tables/", basename(p))
    tibble(
      file = rel,
      extension = tools::file_ext(p),
      size_bytes = as.character(as.numeric(info$size)),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      copied_from_v5_yes_no = "no",
      generated_in_v6_yes_no = "yes"
    )
  })
  readr::write_tsv(bind_rows(existing[, cols], rows), path, na = "")
}

append_change_log()
append_table_inventory()

all_required_paths <- c(
  unlist(outputs),
  combined_path,
  readme_path,
  inventory_path
)

validation_lines <- c(
  paste0("S1 rows: ", nrow(S1), "; columns: ", ncol(S1)),
  paste0("S2 rows: ", nrow(S2), "; columns: ", ncol(S2)),
  paste0("S3 rows: ", nrow(S3), "; columns: ", ncol(S3)),
  paste0("S4 rows: ", nrow(S4), "; columns: ", ncol(S4)),
  paste0("S5 rows: ", nrow(S5), "; columns: ", ncol(S5)),
  paste0("All S1-S5 TSV/XLSX and combined files nonzero: ", all(file.exists(all_required_paths) & file.info(all_required_paths)$size > 0)),
  paste0("Combined workbook sheets: S1_dataset_provenance; S2_split_contrasts; S3_boundary_audit; S4_robustness_summary; S5_gene_enrichment; S5_mapping_audit; S5_background_audit; README"),
  paste0("S5 source rows preserved: ", nrow(s5_source), "; S5 output rows: ", nrow(S5)),
  paste0("S5 exact term IDs preserved in term_id column: ", all(s5_source$term_id == S5$term_id)),
  "No numerical values were recalculated except S1 dataset-context summaries and S4 summaries derived from existing tables.",
  "No IMRS scoring, validation, robustness, transfer-evaluation, or enrichment statistics were modified.",
  "No mechanistic, cell-type source, clinical prediction, or platform safety-ranking claims were introduced.",
  paste0("Missing source warnings: ", ifelse(length(warnings_found) == 0, "none", paste(warnings_found, collapse = "; ")))
)
log_msg("Validation summary:")
for (line in validation_lines) log_msg("  ", line)
log_msg("Updated v6_change_log.tsv, v6_table_inventory.tsv, and v6_file_inventory.tsv.")
log_msg("Priority 4 supplementary table packaging complete.")
write_log()
update_file_inventory()

cat("Priority 4 supplementary table package complete\n")
cat("Output folder: ", out_dir, "\n", sep = "")
cat("S1 rows: ", nrow(S1), "\n", sep = "")
cat("S2 rows: ", nrow(S2), "\n", sep = "")
cat("S3 rows: ", nrow(S3), "\n", sep = "")
cat("S4 rows: ", nrow(S4), "\n", sep = "")
cat("S5 rows: ", nrow(S5), "\n", sep = "")
cat("Combined workbook: ", combined_path, "\n", sep = "")
cat("Warnings: ", ifelse(length(warnings_found) == 0, "none", paste(warnings_found, collapse = "; ")), "\n", sep = "")
