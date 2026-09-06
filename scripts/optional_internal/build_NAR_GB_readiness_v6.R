#!/usr/bin/env Rscript
# Priority 6 NAR Genomics and Bioinformatics readiness package for IMRS v6.
# This script documents and inventories existing inputs/outputs only. It does
# not rerun scoring, validation, robustness, enrichment, figure generation, or
# supplementary table generation.

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

this_file <- imrs_detect_script_path("build_NAR_GB_readiness_v6.R")
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
v6_root <- imrs_config_field_path(config, "manuscript_output_dir")
out_dir <- imrs_config_field_path(config, "internal_readiness_dir", "results_release_templates/internal_qc/NAR_GB_readiness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

delta <- intToUtf8(0x0394)
metric_label <- paste0("delivery-minus-control ", delta, "IMRSz")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
log_path <- file.path(out_dir, "Priority6_NAR_GB_readiness_log.txt")
log_lines <- character()
warnings_found <- character()

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}
add_warning <- function(...) {
  msg <- paste(..., collapse = "")
  warnings_found <<- c(warnings_found, msg)
  log_msg("WARNING: ", msg)
}
write_log <- function() {
  writeLines(log_lines, log_path, useBytes = TRUE)
}

required_packages <- c("readr", "dplyr", "tidyr", "stringr", "tibble", "purrr")
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
})

log_msg("Starting Priority 6 NAR G&B readiness packaging.")
log_msg("Project root: ", project_root)
log_msg("Output folder: ", out_dir)
log_msg("Documentation-only run; no analysis values are modified.")

rel <- function(path) {
  stringr::str_replace(normalizePath(path, winslash = "/", mustWork = FALSE),
                       paste0("^", stringr::fixed(project_root), "/?"), "")
}

read_tsv_chr <- function(path) {
  readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
}

write_tsv <- function(tbl, filename) {
  path <- file.path(out_dir, filename)
  readr::write_tsv(tbl, path, na = "")
  path
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (length(x) == 0) return("NA")
  paste(sort(unique(x)), collapse = "; ")
}

as_num <- function(x) suppressWarnings(as.numeric(x))

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

paths <- list(
  provenance = imrs_config_field_path(config, "provenance_table"),
  role_table = imrs_config_field_path(config, "manuscript_role_table"),
  boundary_audit = imrs_config_field_path(config, "boundary_audit_table"),
  s1 = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "Supplementary_Table_S1_dataset_level_provenance.tsv"),
  s2 = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "Supplementary_Table_S2_split_level_contrasts.tsv"),
  s3 = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "Supplementary_Table_S3_late_context_shifted_boundary_audit.tsv"),
  s4 = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "Supplementary_Table_S4_robustness_summary.tsv"),
  s5 = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "Supplementary_Table_S5_gene_enrichment_results.tsv"),
  supplementary_table_inventory = file.path(imrs_config_field_path(config, "supplementary_tables_dir"), "supplementary_table_inventory.tsv"),
  output_manifest = file.path(v6_root, "v6_output_manifest.tsv"),
  file_inventory = file.path(v6_root, "v6_file_inventory.tsv"),
  table_inventory = file.path(v6_root, "v6_table_inventory.tsv"),
  script_inventory = file.path(v6_root, "v6_script_inventory.tsv"),
  change_log = file.path(v6_root, "v6_change_log.tsv"),
  gene_weights = imrs_config_field_path(config, "frozen_gene_weights"),
  gene_symbol_mapping = imrs_config_field_path(config, "gene_symbol_mapping"),
  step09_eval = imrs_config_field_path(config, "step09_split_eval", "data_release_templates/derived/step09_split_eval.tsv"),
  step09_summary = imrs_config_field_path(config, "step09_split_summary", "data_release_templates/derived/step09_split_summary.tsv"),
  priority3_script = file.path(imrs_config_field_path(config, "active_scripts_dir", "scripts/active_manuscript"), "run_IMRS_retained_gene_enrichment_v6.R"),
  priority4_script = file.path(imrs_config_field_path(config, "active_scripts_dir", "scripts/active_manuscript"), "build_supplementary_tables_v6.R")
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) add_warning("Expected source not found: ", nm, " -> ", paths[[nm]])
}

provenance <- if (file.exists(paths$provenance)) read_tsv_chr(paths$provenance) else tibble()
s1 <- if (file.exists(paths$s1)) read_tsv_chr(paths$s1) else tibble()
s2 <- if (file.exists(paths$s2)) read_tsv_chr(paths$s2) else tibble()
s3 <- if (file.exists(paths$s3)) read_tsv_chr(paths$s3) else tibble()
s4 <- if (file.exists(paths$s4)) read_tsv_chr(paths$s4) else tibble()
s5 <- if (file.exists(paths$s5)) read_tsv_chr(paths$s5) else tibble()
v6_files <- if (file.exists(paths$file_inventory)) read_tsv_chr(paths$file_inventory) else tibble()
v6_tables <- if (file.exists(paths$table_inventory)) read_tsv_chr(paths$table_inventory) else tibble()
v6_scripts <- if (file.exists(paths$script_inventory)) read_tsv_chr(paths$script_inventory) else tibble()
v6_outputs <- if (file.exists(paths$output_manifest)) read_tsv_chr(paths$output_manifest) else tibble()

log_msg("Loaded provenance rows: ", nrow(provenance))
log_msg("Loaded S1-S5 rows: ", paste(c(nrow(s1), nrow(s2), nrow(s3), nrow(s4), nrow(s5)), collapse = ", "))

accession_table <- provenance %>%
  mutate(
    manuscript_analysis_group = controlled_group(.data$manuscript_interpretation_group),
    included = ifelse(toupper(.data$used_for_primary_claim) == "TRUE", "yes", "no"),
    public_repository = dplyr::case_when(
      !is.na(.data$source_database) & nzchar(.data$source_database) & .data$source_database != "NA" ~ .data$source_database,
      stringr::str_detect(.data$dataset_id, "^GSE") | stringr::str_detect(.data$gse_id, "^GSE") ~ "GEO",
      TRUE ~ "NA"
    )
  ) %>%
  group_by(.data$dataset_id, .data$manuscript_analysis_group) %>%
  summarise(
    GEO_accession = collapse_unique(.data$accession),
    PMID = collapse_unique(.data$PMID),
    DOI = collapse_unique(.data$DOI),
    publication_title = collapse_unique(.data$paper_title),
    public_repository = collapse_unique(.data$public_repository),
    data_type = "bulk RNA-seq/count matrix with verified metadata and split-design definitions",
    species = collapse_unique(.data$organism),
    tissue_or_sample_context = collapse_unique(.data$tissue),
    timepoint_h_min = ifelse(all(is.na(as_num(.data$time_h))), "NA", as.character(min(as_num(.data$time_h), na.rm = TRUE))),
    timepoint_h_max = ifelse(all(is.na(as_num(.data$time_h))), "NA", as.character(max(as_num(.data$time_h), na.rm = TRUE))),
    delivery_platform_or_modality = collapse_unique(.data$delivery_platform),
    role_in_IMRS = collapse_unique(.data$pipeline_role),
    included_in_primary_claim_yes_no = ifelse(any(.data$included == "yes", na.rm = TRUE), "yes", "no"),
    raw_data_access_note = "Raw data are public through the listed accession when available; derived metadata and split definitions are documented in the release package.",
    derived_data_available_in_package_yes_no = "yes",
    source_file = rel(paths$provenance),
    notes = dplyr::case_when(
      PMID == "NA" | DOI == "NA" ~ "PMID and/or DOI unavailable in local provenance source.",
      TRUE ~ "Public dataset provenance captured from existing local audit source."
    ),
    .groups = "drop"
  ) %>%
  transmute(
    dataset_id, GEO_accession, PMID, DOI, publication_title, public_repository,
    data_type, species, tissue_or_sample_context, timepoint_h_min, timepoint_h_max,
    delivery_platform_or_modality, manuscript_analysis_group, role_in_IMRS,
    included_in_primary_claim_yes_no, raw_data_access_note,
    derived_data_available_in_package_yes_no, source_file, notes
  ) %>%
  arrange(.data$manuscript_analysis_group, .data$dataset_id)
accession_path <- write_tsv(accession_table, "IMRS_public_dataset_accession_table_v6.tsv")

code_statement <- paste(
  "Analysis code for the frozen, transfer-oriented IMRS scoring framework, figure generation, robustness analyses, comparator-signature benchmarking, retained-gene enrichment analysis, and supplementary table generation is available at https://github.com/YanshengLuo/InnateImmuneResponseScore. The repository includes scripts used for metadata processing, locked-anchor coefficient construction, frozen scoring, delivery-versus-control split evaluation, manuscript role curation, robustness checks, figure generation, and supplementary table packaging. Versioned release information and software/session details are provided in the reproducibility materials.",
  "",
  "For peer review, the repository release candidate and accompanying derived input/output package correspond to the v6 manuscript-output version. No archival DOI has been assigned in this release candidate.",
  sep = "\n"
)
writeLines(code_statement, file.path(out_dir, "IMRS_code_availability_statement_candidate.txt"), useBytes = TRUE)

data_statement <- paste(
  paste0("All raw sequencing datasets analyzed in this study were obtained from public repositories, primarily GEO, under the accession numbers listed in Supplementary Table S1. Reprocessed metadata, split-design files, frozen IMRS coefficient tables, ", metric_label, " summaries, robustness outputs, enrichment outputs, and manuscript-ready supplementary tables are available in the repository release candidate at https://github.com/YanshengLuo/InnateImmuneResponseScore. The derived data package includes dataset-level provenance, split-level contrast definitions, manuscript role assignments, late or context-shifted boundary-context annotations, robustness summaries, and retained-gene enrichment results. No new human-subject data were generated for this study."),
  "",
  "Raw sequencing files should be retrieved from GEO or other indicated public repositories by accession rather than redistributed in this source repository unless redistribution is separately confirmed as appropriate.",
  sep = "\n"
)
writeLines(data_statement, file.path(out_dir, "IMRS_data_availability_statement_candidate.txt"), useBytes = TRUE)

known_main_outputs <- tibble(
  manuscript_component = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Supplementary Figure S1", "Supplementary Figure S2"),
  figure_or_table_generated = c(
    "Figure1_main_v5.png/pdf/svg; Figure1A_IMRS_merged_workflow_v5.png/pdf/svg; Figure1B_dataset_tissue_response_landscape_v5_corrected.png/pdf/svg",
    "Figure2_main_v5.png/pdf/svg and intermediate Figure2A-D panel outputs",
    "Figure3_main_v5.png/pdf/svg and intermediate Figure3A-B panel outputs",
    "Figure4_main_v5.png/pdf/svg and intermediate Figure4A-B panel outputs",
    "Figure5_main_v5.png/pdf/svg and intermediate Figure5A-D panel outputs",
    "FigureS_comparator_benchmarking_v5.png/pdf/svg; FigureS_validation_faceted_summary_v5.png/pdf/svg; FigureS_weak_context_interpretation_categories_v5.png/pdf/svg",
    "Priority3_gene_program_enrichment/figures/FigureS2A-C and combined outputs"
  ),
  output_files_detected = c(
    "Figure1_main_v5.*; Figure1A_IMRS_merged_workflow_v5.*; Figure1B_dataset_tissue_response_landscape_v5_corrected.*",
    "Figure2_main_v5.*; intermediate_panels/Figure2A-D_v5.*",
    "Figure3_main_v5.*; intermediate_panels/Figure3A-B_v5.*",
    "Figure4_main_v5.*; intermediate_panels/Figure4A-B_v5.*",
    "Figure5_main_v5.*; intermediate_panels/Figure5A-D_v5.*",
    "FigureS_comparator_benchmarking_v5.*; FigureS_validation_faceted_summary_v5.*; FigureS_weak_context_interpretation_categories_v5.*",
    "Priority3_gene_program_enrichment/figures/FigureS2A-C_*; FigureS2_gene_program_enrichment_combined.*"
  ),
  script_file = c(
    "scripts/R/build_merged_imrs_workflow_v5.R; scripts/01_regenerate_changed_v6.R",
    "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
    "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
    "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
    "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
    "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
    "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R"
  ),
  script_purpose = c(
    "Build and regenerate Figure 1 workflow/landscape outputs.",
    "Generate Figure 2 retained-gene/coefficient/support panels.",
    "Generate Figure 3 transfer validation panels.",
    "Generate Figure 4 boundary-context panels.",
    "Generate Figure 5 robustness/non-degeneracy panels.",
    "Generate comparator, validation, and context-shifted supplementary panels.",
    "Generate retained-gene program enrichment tables and Figure S2."
  ),
  input_files_detected = c(
    "v6 scripts/tables plus Step09 and provenance outputs",
    "05_score/anchors/gene_weights.tsv; 05_score/anchors/support_by_dataset.tsv; audit role tables",
    "05_score/transfer/eval/step09_split_eval.tsv; audit role/provenance tables",
    "audit/results/weak_dataset_paper_context_audit.tsv; Step09 outputs",
    "05_score/publication_extra_analyses/results; audit/results/leave_one_gene_out_summary.tsv; gene_dominance_summary.tsv",
    "05_score/publication_extra_analyses/results/baseline_signature_*; audit validation/context outputs",
    "05_score/anchors/gene_weights.tsv; gene_symbol_mapping.tsv; gene_power.tsv"
  ),
  runnable_status_inferred = c(
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable inferred from existing v6 generation logs, but not rerun by Priority 6",
    "runnable checked in prior Priority 3 runs; not rerun by Priority 6"
  ),
  last_modified_time = "see file inventory",
  notes = "Mapping is reviewer-facing documentation only; no analysis values were modified."
)

local_scripts <- list.files(imrs_config_field_path(config, "active_scripts_dir", "scripts/active_manuscript"),
                            pattern = "\\.(R|ps1)$", recursive = TRUE, full.names = TRUE)
pipeline_scripts <- c(
  list.files(imrs_config_field_path(config, "audit_scripts_dir", "scripts/full_pipeline"),
             pattern = "\\.R$", recursive = TRUE, full.names = TRUE),
  list.files(imrs_config_field_path(config, "publication_extra_scripts_dir", "scripts/full_pipeline"),
             pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
)
script_candidates <- unique(c(local_scripts, pipeline_scripts))
script_rows <- tibble(
  script_file = rel(script_candidates),
  script_purpose = dplyr::case_when(
    str_detect(script_file, "Priority3_gene_program_enrichment") ~ "Retained-gene enrichment analysis and Figure S2 generation.",
    str_detect(script_file, "supplementary_tables") ~ "Supplementary table packaging.",
    str_detect(script_file, "NAR_GB_readiness") ~ "NAR G&B readiness documentation packaging.",
    str_detect(script_file, "00_generate_all_reorganized_figures") ~ "Main and supplementary v5/v6 figure generation.",
    str_detect(script_file, "01_regenerate_changed_v6") ~ "Targeted v6 figure regeneration.",
    str_detect(script_file, "02_create_v6_metadata") ~ "v6 metadata/inventory creation.",
    str_detect(script_file, "label_permutation") ~ "Label permutation robustness analysis.",
    str_detect(script_file, "baseline_signature") ~ "Comparator-signature benchmarking.",
    str_detect(script_file, "coefficient_sensitivity") ~ "Coefficient sensitivity summary.",
    str_detect(script_file, "leave_one_gene") ~ "Leave-one-gene-out robustness analysis.",
    str_detect(script_file, "gene_dominance") ~ "Gene dominance robustness analysis.",
    str_detect(script_file, "threshold_sensitivity") ~ "Threshold sensitivity analysis.",
    str_detect(script_file, "leave_one_anchor") ~ "Leave-one-anchor-out robustness analysis.",
    TRUE ~ "Existing project/audit script; purpose inferred from filename."
  ),
  input_files_detected = "See script body and source manifests; not rerun in Priority 6.",
  output_files_detected = "Existing outputs documented in v6 manifests and readiness package.",
  figure_or_table_generated = dplyr::case_when(
    str_detect(script_file, "Priority3_gene_program_enrichment") ~ "Supplementary Figure S2; Supplementary Table S5 enrichment source outputs",
    str_detect(script_file, "supplementary_tables") ~ "Supplementary Tables S1-S5 package",
    str_detect(script_file, "NAR_GB_readiness") ~ "NAR G&B readiness package",
    str_detect(script_file, "generate_all|regenerate_changed|build_merged") ~ "Main/supplementary figures",
    TRUE ~ "Audit/robustness/provenance outputs"
  ),
  manuscript_component = dplyr::case_when(
    str_detect(script_file, "Priority3_gene_program_enrichment") ~ "Supplementary Figure S2 and gene-program enrichment results",
    str_detect(script_file, "supplementary_tables") ~ "Supplementary Tables S1-S5",
    str_detect(script_file, "NAR_GB_readiness") ~ "Reproducibility and availability materials",
    str_detect(script_file, "publication_extra|audit/scripts") ~ "Robustness/provenance/context analyses",
    TRUE ~ "Main and supplementary figure generation"
  ),
  last_modified_time = format(file.info(script_candidates)$mtime, "%Y-%m-%d %H:%M:%S"),
  runnable_status_inferred = "not run by Priority 6; status inferred from existing outputs or filename",
  notes = "Script included for reviewer-facing traceability."
)
script_manifest <- bind_rows(known_main_outputs %>% select(names(script_rows)), script_rows) %>%
  distinct(script_file, figure_or_table_generated, .keep_all = TRUE) %>%
  arrange(.data$manuscript_component, .data$script_file)
script_manifest_path <- write_tsv(script_manifest, "IMRS_script_to_output_manifest_v6.tsv")

input_candidates <- c(
  list.files(imrs_config_field_path(config, "curated_metadata_dir"), pattern = "\\.(tsv|csv|txt)$", recursive = TRUE, full.names = TRUE),
  list.files(imrs_config_field_path(config, "derived_data_dir"), pattern = "\\.(tsv|csv)$", recursive = TRUE, full.names = TRUE),
  list.files(imrs_config_field_path(config, "priority3_enrichment_dir"), pattern = "\\.(tsv|xlsx)$", recursive = TRUE, full.names = TRUE),
  list.files(imrs_config_field_path(config, "supplementary_tables_dir"), pattern = "\\.(tsv|xlsx|txt)$", recursive = TRUE, full.names = TRUE)
)
input_candidates <- unique(input_candidates[file.exists(input_candidates)])
input_manifest <- tibble(input_file_abs = input_candidates) %>%
  mutate(
    input_file = rel(.data$input_file_abs),
    file_type = tools::file_ext(.data$input_file_abs),
    inferred_purpose = case_when(
      str_detect(input_file, "verified_metadata/.*/scoring") ~ "verified metadata file for scoring",
      str_detect(input_file, "verified_metadata/.*/splited|verified_metadata/splited") ~ "split-design file for delivery-versus-control contrast",
      str_detect(input_file, "gene_weights") ~ "frozen IMRS coefficient/gene-weight table",
      str_detect(input_file, "gene_symbol_mapping") ~ "gene identifier mapping table",
      str_detect(input_file, "gene_power|gene_heterogeneity|support_by_dataset") ~ "anchor construction support/background table",
      str_detect(input_file, "step09_split") ~ "split-level transfer-evaluation output used as source input",
      str_detect(input_file, "supplement_dataset_split_provenance") ~ "provenance and split-context source",
      str_detect(input_file, "weak_dataset|boundary") ~ "late or context-shifted boundary-context audit source",
      str_detect(input_file, "Supplementary_Table_S5_IMRS_gene_enrichment") ~ "Priority 3 retained-gene enrichment source",
      str_detect(input_file, "IMRS_retained_gene_mapping|IMRS_enrichment_background") ~ "Priority 3 mapping/background audit",
      TRUE ~ "derived source/support table"
    ),
    dataset_id_if_applicable = str_extract(input_file, "GSE[0-9]+(?:_[A-Za-z0-9_]+)?"),
    dataset_id_if_applicable = ifelse(is.na(dataset_id_if_applicable), "NA", dataset_id_if_applicable),
    public_or_derived = case_when(
      str_detect(input_file, "00_metadata|raw|counts|featurecounts") ~ "public-derived/reprocessed",
      str_detect(input_file, "05_score|audit|final_manuscripts") ~ "derived",
      TRUE ~ "unknown"
    ),
    source_repository_or_project_path = input_file,
    expected_in_release_yes_no = ifelse(str_detect(input_file, "raw|counts") & !str_detect(input_file, "verified_metadata"), "no", "yes"),
    redistribution_allowed_unknown_yes_no = ifelse(str_detect(input_file, "raw|counts|featurecounts"), "yes", "no"),
    notes = case_when(
      str_detect(input_file, "raw|counts|featurecounts") ~ "Raw or reprocessed count data may need accession-based retrieval rather than redistribution.",
      str_detect(input_file, "verified_metadata") ~ "Metadata/split design should be released for reproducibility where permitted.",
      TRUE ~ "Derived project input/support file for reviewer traceability."
    )
  ) %>%
  select(-input_file_abs) %>%
  arrange(.data$input_file)
input_manifest_path <- write_tsv(input_manifest, "IMRS_input_data_manifest_v6.tsv")

derived_files <- list.files(v6_root, recursive = TRUE, full.names = TRUE)
derived_files <- derived_files[file.exists(derived_files) & !dir.exists(derived_files)]
derived_manifest <- tibble(output_abs = derived_files) %>%
  mutate(
    output_file = rel(.data$output_abs),
    file_type = tools::file_ext(.data$output_abs),
    figure_or_table_id = case_when(
      str_detect(output_file, "Figure1|Figure1A|Figure1B") ~ "Figure 1",
      str_detect(output_file, "Figure2") ~ "Figure 2",
      str_detect(output_file, "Figure3") ~ "Figure 3",
      str_detect(output_file, "Figure4") ~ "Figure 4",
      str_detect(output_file, "Figure5") ~ "Figure 5",
      str_detect(output_file, "FigureS2|Priority3_gene_program_enrichment") ~ "Supplementary Figure S2 / Priority 3",
      str_detect(output_file, "FigureS_|Supplementary_FigureS1") ~ "Supplementary Figure S1",
      str_detect(output_file, "Supplementary_Table_S1") ~ "Supplementary Table S1",
      str_detect(output_file, "Supplementary_Table_S2") ~ "Supplementary Table S2",
      str_detect(output_file, "Supplementary_Table_S3") ~ "Supplementary Table S3",
      str_detect(output_file, "Supplementary_Table_S4") ~ "Supplementary Table S4",
      str_detect(output_file, "Supplementary_Table_S5") ~ "Supplementary Table S5",
      str_detect(output_file, "NAR_GB_readiness") ~ "NAR G&B readiness",
      TRUE ~ "manifest/log/supporting output"
    ),
    inferred_purpose = case_when(
      file_type %in% c("png", "pdf", "svg") ~ "figure output",
      str_detect(output_file, "supplementary_tables") ~ "manuscript-ready supplementary table package",
      str_detect(output_file, "Priority3_gene_program_enrichment") ~ "retained-gene enrichment output",
      str_detect(output_file, "manifest|inventory") ~ "manifest/inventory output",
      str_detect(output_file, "log") ~ "generation or validation log",
      TRUE ~ "derived manuscript support output"
    ),
    generated_by_script = case_when(
      str_detect(output_file, "Priority3_gene_program_enrichment") ~ "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R",
      str_detect(output_file, "supplementary_tables") ~ "supplementary_tables/build_supplementary_tables_v6.R",
      str_detect(output_file, "NAR_GB_readiness") ~ "NAR_GB_readiness/build_NAR_GB_readiness_v6.R",
      str_detect(output_file, "Figure1A|Figure1B|Figure[1-5]|FigureS_") ~ "scripts/00_generate_all_reorganized_figures_revised_v5.R; scripts/01_regenerate_changed_v6.R",
      TRUE ~ "see v6 manifests/logs"
    ),
    manuscript_or_supplementary = case_when(
      str_detect(figure_or_table_id, "^Figure [1-5]$") ~ "manuscript",
      str_detect(figure_or_table_id, "Supplementary|NAR") ~ "supplementary/reviewer-facing",
      TRUE ~ "supporting"
    ),
    expected_in_release_yes_no = "yes",
    notes = "Derived output; included for reviewer traceability where useful."
  ) %>%
  select(-output_abs) %>%
  arrange(.data$figure_or_table_id, .data$output_file)
derived_manifest_path <- write_tsv(derived_manifest, "IMRS_derived_output_manifest_v6.tsv")

checklist <- tibble::tribble(
  ~check_item, ~status, ~evidence_file, ~notes,
  "NAR_GB_readiness folder exists", "complete", "NAR_GB_readiness", "Created in v6 folder.",
  "Public dataset accession table", "complete", "NAR_GB_readiness/IMRS_public_dataset_accession_table_v6.tsv", "Generated from existing provenance source.",
  "Code availability statement candidate", "complete", "NAR_GB_readiness/IMRS_code_availability_statement_candidate.txt", "Uses the public repository URL and records that no archival DOI is assigned.",
  "Data availability statement candidate", "complete", "NAR_GB_readiness/IMRS_data_availability_statement_candidate.txt", "Uses the public repository URL and states the raw-data retrieval policy.",
  "Reproducibility README", "complete", "NAR_GB_readiness/IMRS_reproducibility_README_v6.md", "Reviewer-facing workflow and limitation summary.",
  "Script-to-output manifest", "complete", "NAR_GB_readiness/IMRS_script_to_output_manifest_v6.tsv", "Maps major scripts to outputs and manuscript components.",
  "Input data manifest", "complete", "NAR_GB_readiness/IMRS_input_data_manifest_v6.tsv", "Classifies metadata, split designs, coefficients, and derived sources.",
  "Derived output manifest", "complete", "NAR_GB_readiness/IMRS_derived_output_manifest_v6.tsv", "Includes figures, supplementary tables, enrichment outputs, manifests, and logs.",
  "Software/session info", "complete", "NAR_GB_readiness/IMRS_software_session_info_v6.txt", "Captured from this R environment.",
  "Repository release plan", "complete", "NAR_GB_readiness/IMRS_repository_release_plan_v6.md", "Includes GitHub/Zenodo/reviewer-access plan.",
  "No analysis value changes", "complete", "Priority6_NAR_GB_readiness_log.txt", "Documentation-only run."
)
checklist_path <- write_tsv(checklist, "IMRS_reproducibility_checklist_v6.tsv")

repo_plan <- paste(
  "# IMRS Repository Release Plan v6",
  "",
  "## Recommended Repository Structure",
  "",
  "```text",
  "IMRS/",
  "  README.md",
  "  LICENSE",
  "  CITATION.cff",
  "  renv.lock or environment.yml",
  "  data/",
  "    external_accessions/",
  "    metadata/",
  "    derived/",
  "  scripts/",
  "    01_metadata/",
  "    02_design/",
  "    03_anchor_model/",
  "    04_scoring/",
  "    05_transfer_evaluation/",
  "    06_robustness/",
  "    07_figures/",
  "    08_supplementary_tables/",
  "  results/",
  "    figures/",
  "    supplementary_figures/",
  "    supplementary_tables/",
  "    logs/",
  "  docs/",
  "    reproducibility/",
  "    data_dictionary/",
  "```",
  "",
  "## Public GitHub Contents",
  "",
  "Include analysis scripts, metadata-processing scripts, split-design builders, locked-anchor coefficient construction code, frozen scoring code, transfer-evaluation code, robustness scripts, comparator-signature benchmarking, retained-gene enrichment, figure generation, supplementary table packaging, manifests, data dictionaries, and documentation.",
  "",
  "## Zenodo or Equivalent Archive",
  "",
  "Archive a versioned release matching the submitted manuscript. Include derived tables, frozen gene weights, verified metadata and split-design definitions where redistribution is permitted, figure outputs, supplementary tables, reproducibility manifests, logs, and session information.",
  "",
  "## Data Not to Redistribute Directly",
  "",
  "Do not redistribute raw public GEO data if the preferred journal/repository practice is to retrieve them from their original public accessions. Instead, provide accession numbers, retrieval scripts, metadata curation, split definitions, and checksums or provenance notes sufficient to rebuild the analysis inputs.",
  "",
  "## Review Access",
  "",
  "If the repository remains private before acceptance, provide a reviewer-access link or archive. The reviewer package should include the v6 outputs, scripts, manifests, data/code availability statements, and this NAR G&B readiness folder.",
  "",
  "## Scope Statement",
  "",
  paste0("The release documents the frozen, transfer-oriented IMRS scoring framework and ", metric_label, " analyses. It does not frame IMRS as a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool."),
  sep = "\n"
)
writeLines(repo_plan, file.path(out_dir, "IMRS_repository_release_plan_v6.md"), useBytes = TRUE)

review_inventory <- tibble::tribble(
  ~package_item, ~file_or_folder, ~purpose, ~required_for_reproduction_yes_no, ~manuscript_component_supported, ~ready_for_reviewer_yes_no, ~missing_or_needs_cleanup, ~notes,
  "Reproducibility README", "NAR_GB_readiness/IMRS_reproducibility_README_v6.md", "Reviewer-facing overview and workflow.", "yes", "all", "yes", "none", "Generated by Priority 6.",
  "Code availability statement", "NAR_GB_readiness/IMRS_code_availability_statement_candidate.txt", "Manuscript code-availability candidate.", "no", "availability statements", "yes", "none", "Add an archival DOI later only if one is issued.",
  "Data availability statement", "NAR_GB_readiness/IMRS_data_availability_statement_candidate.txt", "Manuscript data-availability candidate.", "no", "availability statements", "yes", "none", "Raw GEO data are referenced by accession rather than redistributed.",
  "Public dataset accessions", "NAR_GB_readiness/IMRS_public_dataset_accession_table_v6.tsv", "Public dataset provenance and roles.", "yes", "Supplementary Table S1/data availability", "yes", "PMID/DOI unavailable for rows where local provenance lacks them", "Do not invent missing values.",
  "Script-to-output manifest", "NAR_GB_readiness/IMRS_script_to_output_manifest_v6.tsv", "Map scripts to manuscript outputs.", "yes", "figures/tables/robustness/enrichment", "yes", "runnable status is inferred unless previously checked", "Reviewer traceability.",
  "Input data manifest", "NAR_GB_readiness/IMRS_input_data_manifest_v6.tsv", "Classify inputs and release expectations.", "yes", "data provenance", "yes", "raw-data redistribution status may require repository policy decision", "Accession-based retrieval expected for public raw data.",
  "Derived output manifest", "NAR_GB_readiness/IMRS_derived_output_manifest_v6.tsv", "Inventory derived outputs.", "yes", "all", "yes", "none", "Includes figures/tables/logs/manifests.",
  "Software session info", "NAR_GB_readiness/IMRS_software_session_info_v6.txt", "R/session/package versions.", "yes", "reproducibility", "yes", "none", "Generated from local R session.",
  "Repository release plan", "NAR_GB_readiness/IMRS_repository_release_plan_v6.md", "Plan GitHub/Zenodo/reviewer release.", "no", "reproducibility/data/code availability", "yes", "requires repository setup", "Manual action before submission.",
  "Main figures", "Figure1-Figure5 PNG/PDF/SVG outputs", "Manuscript figure outputs.", "yes", "Figures 1-5", "yes", "none", "Existing v6 outputs inventoried.",
  "Supplementary figures", "Supplementary Figure S1/S2 outputs", "Supplementary figure outputs.", "yes", "Supplementary Figures S1-S2", "yes", "none", "Existing v6 outputs inventoried.",
  "Supplementary tables", "supplementary_tables/", "Manuscript-ready S1-S5 package.", "yes", "Supplementary Tables S1-S5", "yes", "none", "Generated by Priority 4.",
  "Logs", "v6 logs and Priority6_NAR_GB_readiness_log.txt", "Generation and validation trace.", "yes", "reviewer-facing reproducibility", "yes", "none", "Documentation-only log generated."
)
review_inventory_path <- write_tsv(review_inventory, "IMRS_review_package_inventory_v6.tsv")

gap_audit <- tibble::tribble(
  ~requirement_or_expectation, ~current_status, ~evidence_file, ~gap_level, ~recommended_action, ~notes,
  "Public dataset accession table", "Generated from existing provenance.", "NAR_GB_readiness/IMRS_public_dataset_accession_table_v6.tsv", "complete", "Review NA PMID/DOI rows before submission.", "Do not invent unavailable values.",
  "Code availability statement", "Candidate includes public repository URL.", "NAR_GB_readiness/IMRS_code_availability_statement_candidate.txt", "complete", "Add an archival DOI only if one is issued.", "No DOI asserted in this release candidate.",
  "Data availability statement", "Candidate includes repository URL and raw-data access policy.", "NAR_GB_readiness/IMRS_data_availability_statement_candidate.txt", "complete", "Keep accession table synchronized with the manuscript.", "Raw GEO data should be accessed by accession when redistribution is not appropriate.",
  "Reproducibility README", "Created.", "NAR_GB_readiness/IMRS_reproducibility_README_v6.md", "complete", "Review against final repository structure.", "Describes workflow and scope.",
  "Script-to-output manifest", "Created.", "NAR_GB_readiness/IMRS_script_to_output_manifest_v6.tsv", "partial", "Confirm runnable status during repository release dry run.", "Priority 6 did not rerun scripts.",
  "Input data manifest", "Created.", "NAR_GB_readiness/IMRS_input_data_manifest_v6.tsv", "partial", "Decide which raw/reprocessed count files can be redistributed.", "Accessions and rebuild scripts should cover non-redistributed public data.",
  "Derived output manifest", "Created.", "NAR_GB_readiness/IMRS_derived_output_manifest_v6.tsv", "complete", "Keep synchronized with final v6 release.", "Includes figures, tables, manifests, and logs.",
  "Software/session info", "Captured from local R session.", "NAR_GB_readiness/IMRS_software_session_info_v6.txt", "complete", "Consider freezing with renv.lock or environment.yml for release.", "Major packages listed.",
  "Supplementary tables S1-S5", "Generated in Priority 4.", "supplementary_tables/", "complete", "Include in release archive.", "No values changed in Priority 6.",
  "Figure role manifest", "Existing v6 output and file manifests present.", "v6_output_manifest.tsv; figure_v6_manifest.tsv", "complete", "Review final figure filenames before release.", "Manual captions may still require manuscript assembly.",
  "Controlled terminology consistency", "Priority 6 package uses controlled terms.", "NAR_GB_readiness/IMRS_reproducibility_README_v6.md", "complete", "Run final manuscript-wide terminology audit before submission.", "Avoids unsupported claims.",
  "Clear scope statement", "Included in README and release plan.", "NAR_GB_readiness/IMRS_reproducibility_README_v6.md; IMRS_repository_release_plan_v6.md", "complete", "Keep in final manuscript text.", "States IMRS is not a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool.",
  "Data/code repository release plan", "Created.", "NAR_GB_readiness/IMRS_repository_release_plan_v6.md", "manual_action_needed", "Create repository, license, citation, environment lockfile, and archive DOI.", "Plan is ready, external repository not created by this task.",
  "Reviewer access package", "Inventory created.", "NAR_GB_readiness/IMRS_review_package_inventory_v6.tsv", "manual_action_needed", "Provide private repository link or reviewer archive.", "Needs manual upload/access control."
)
gap_audit_path <- write_tsv(gap_audit, "IMRS_NAR_GB_submission_gap_audit_v6.tsv")

pkg_versions <- c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "patchwork", "scales",
                  "openxlsx", "writexl", "clusterProfiler", "org.Mm.eg.db", "ReactomePA",
                  "msigdbr", "AnnotationDbi", "DESeq2")
session_path <- file.path(out_dir, "IMRS_software_session_info_v6.txt")
session_lines <- capture.output({
  cat("IMRS v6 software/session information\n")
  cat("Generated: ", timestamp, "\n\n", sep = "")
  cat("R sessionInfo():\n")
  print(sessionInfo())
  cat("\nMajor package versions:\n")
  for (pkg in pkg_versions) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(pkg, "\t", as.character(utils::packageVersion(pkg)), "\n", sep = "")
    } else {
      cat(pkg, "\tunavailable in this R library\n", sep = "")
    }
  }
})
writeLines(session_lines, session_path, useBytes = TRUE)

readme <- paste(
  "# IMRS v6 Reproducibility README",
  "",
  "## 1. Project Overview",
  "",
  paste0("IMRS is a frozen, transfer-oriented transcriptomic scoring framework for acute delivery-associated innate transcriptional responses. The v6 manuscript materials document locked-anchor coefficient construction, frozen scoring, split-level delivery-versus-control evaluation, robustness checks, comparator analyses, gene-program enrichment, and manuscript-ready figure/table outputs. IMRS is not a mechanistic pathway model, clinical reactogenicity predictor, or delivery-platform safety ranking tool."),
  "",
  "## 2. Folder Structure",
  "",
  "- `revised_plots_v6/`: manuscript-facing figures, manifests, logs, scripts, and current v6 outputs.",
  "- `revised_plots_v6/intermediate_panels/`: figure panel-level PNG/PDF/SVG outputs.",
  "- `revised_plots_v6/Priority3_gene_program_enrichment/`: retained-gene enrichment script, tables, Figure S2 outputs, captions, and logs.",
  "- `revised_plots_v6/supplementary_tables/`: manuscript-ready Supplementary Tables S1-S5, combined workbook, README, inventory, and generation log.",
  "- `revised_plots_v6/NAR_GB_readiness/`: reproducibility README, availability statements, accession table, manifests, session info, release plan, review inventory, and gap audit.",
  "",
  "## 3. Data Sources",
  "",
  "Public datasets are listed in `IMRS_public_dataset_accession_table_v6.tsv`, including GEO accession, PMID, DOI, publication title, species, tissue/sample context, time range, delivery modality, and manuscript analysis group. Controlled manuscript analysis groups are Locked anchor, Primary acute validation, Extended validation, Secondary support, and Excluded/unclear when local metadata require caution.",
  "",
  "## 4. Reproducibility Workflow",
  "",
  "1. Verified metadata and raw/reprocessed count matrices define delivery-versus-control split designs.",
  "2. Locked-anchor datasets are used for locked-anchor coefficient construction.",
  "3. Frozen coefficients are applied without refitting in the frozen, transfer-oriented IMRS scoring framework.",
  paste0("4. Split-level delivery-versus-control evaluation produces ", metric_label, " and secondary AUC summaries."),
  "5. Manuscript role cleanup assigns Locked anchor, Primary acute validation, Extended validation, Secondary support, or Excluded/unclear labels.",
  "6. Boundary-context audit retains late or context-shifted settings as boundary-setting evidence.",
  "7. Robustness checks evaluate score stability and non-degeneracy, including label permutation, leave-one-gene-out, gene dominance, threshold sensitivity, leave-one-anchor-out, comparator benchmarking, and coefficient sensitivity.",
  "8. Comparator signatures contextualize related immune-signature behavior without replacing IMRS.",
  "9. Gene-program enrichment provides program-level context for the acute delivery-associated innate transcriptional program represented by retained frozen coefficients.",
  "",
  "## 5. Main Figure Generation",
  "",
  "- Figure 1: `Figure1_main_v5.*`, `Figure1A_IMRS_merged_workflow_v5.*`, and `Figure1B_dataset_tissue_response_landscape_v5_corrected.*`; source scripts include `scripts/R/build_merged_imrs_workflow_v5.R` and `scripts/01_regenerate_changed_v6.R`.",
  "- Figure 2: `Figure2_main_v5.*` and `intermediate_panels/Figure2A-D_v5.*`; source script mapping is recorded in `IMRS_script_to_output_manifest_v6.tsv` and v6 figure manifests.",
  "- Figure 3: `Figure3_main_v5.*` and `intermediate_panels/Figure3A-B_v5.*`.",
  "- Figure 4: `Figure4_main_v5.*` and `intermediate_panels/Figure4A-B_v5.*`.",
  "- Figure 5: `Figure5_main_v5.*` and `intermediate_panels/Figure5A-D_v5.*`.",
  "- Supplementary Figure S1: comparator, validation, and context-shifted supplementary figure outputs in the v6 root and intermediate panel folders.",
  "- Supplementary Figure S2: `Priority3_gene_program_enrichment/figures/FigureS2A-C_*` and `FigureS2_gene_program_enrichment_combined.*`, generated by `Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R`.",
  "",
  "## 6. Supplementary Table Generation",
  "",
  "- Supplementary Table S1: dataset-level provenance from `audit/results/supplement_dataset_split_provenance_v7.tsv`.",
  "- Supplementary Table S2: split-level contrast table from the same provenance source.",
  "- Supplementary Table S3: late/context-shifted boundary audit from `audit/results/weak_dataset_paper_context_audit.tsv`.",
  "- Supplementary Table S4: robustness summary from existing robustness outputs in `audit/results/` and `05_score/publication_extra_analyses/results/`.",
  "- Supplementary Table S5: retained-gene enrichment results from `Priority3_gene_program_enrichment/tables/Supplementary_Table_S5_IMRS_gene_enrichment_all.tsv`.",
  "- Combined workbook and inventories are generated by `supplementary_tables/build_supplementary_tables_v6.R`.",
  "",
  "## 7. Software/Package Requirements",
  "",
  "R session information and major package versions are captured in `IMRS_software_session_info_v6.txt`. Major packages include readr, dplyr, tidyr, stringr, ggplot2, patchwork, scales, writexl/openxlsx, clusterProfiler, org.Mm.eg.db, ReactomePA, msigdbr, AnnotationDbi, and DESeq2 if available.",
  "",
  "## 8. Reproduction Instructions",
  "",
  "1. Restore the project folder structure and install the R package environment documented in `IMRS_software_session_info_v6.txt`.",
  "2. Retrieve public raw data from the accessions listed in `IMRS_public_dataset_accession_table_v6.tsv` where raw matrices are not redistributed.",
  "3. Use verified metadata and split-design files recorded in `IMRS_input_data_manifest_v6.tsv`.",
  "4. Run the IMRS pipeline scripts for metadata processing, locked-anchor coefficient construction, frozen scoring, and split-level transfer evaluation if rebuilding from raw inputs.",
  "5. Run robustness and comparator scripts listed in `IMRS_script_to_output_manifest_v6.tsv` if rebuilding robustness outputs.",
  "6. Run v6 figure-generation scripts to regenerate main and supplementary figure outputs.",
  "7. Run `Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R` to regenerate retained-gene enrichment outputs.",
  "8. Run `supplementary_tables/build_supplementary_tables_v6.R` to regenerate Supplementary Tables S1-S5.",
  "9. Run `NAR_GB_readiness/build_NAR_GB_readiness_v6.R` to regenerate these readiness materials.",
  "",
  "Caption note: final manuscript captions and some assembled figure subtitles may be manually controlled in manuscript/PowerPoint/Word sources. Manual caption replacements are tracked in `v6_manual_caption_replacement_table.tsv` when no editable generated source is present.",
  "",
  "## 9. Data/Code Availability",
  "",
  "Candidate statements are provided in `IMRS_code_availability_statement_candidate.txt` and `IMRS_data_availability_statement_candidate.txt`.",
  "",
  "## 10. Known Limitations",
  "",
  "- Several split contrasts have small sample sizes.",
  "- AUC is secondary and is not framed as clinical prediction.",
  "- Current evidence is mostly mouse bulk RNA-seq, with context-specific limitations noted in provenance and boundary-audit tables.",
  "- Late or context-shifted settings are interpreted as boundary-setting evidence rather than primary validation.",
  "- No clinical reactogenicity endpoint is modeled.",
  "- Gene-program enrichment provides program-level biological context and does not establish causal pathways or cell-type sources.",
  sep = "\n"
)
writeLines(readme, file.path(out_dir, "IMRS_reproducibility_README_v6.md"), useBytes = TRUE)

update_change_log <- function() {
  path <- paths$change_log
  cols <- c("file", "panel_or_figure", "change_type", "old_text", "new_text", "plot_visible_or_caption", "regenerated_or_copied", "notes")
  existing <- if (file.exists(path)) read_tsv_chr(path) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% filter(.data$change_type != "Priority6_NAR_GB_readiness")
  created_files <- c(
    "IMRS_reproducibility_README_v6.md",
    "IMRS_code_availability_statement_candidate.txt",
    "IMRS_data_availability_statement_candidate.txt",
    "IMRS_reproducibility_checklist_v6.tsv",
    "IMRS_public_dataset_accession_table_v6.tsv",
    "IMRS_script_to_output_manifest_v6.tsv",
    "IMRS_input_data_manifest_v6.tsv",
    "IMRS_derived_output_manifest_v6.tsv",
    "IMRS_software_session_info_v6.txt",
    "IMRS_repository_release_plan_v6.md",
    "IMRS_review_package_inventory_v6.tsv",
    "IMRS_NAR_GB_submission_gap_audit_v6.tsv",
    "Priority6_NAR_GB_readiness_log.txt"
  )
  rows <- tibble(
    file = paste0("NAR_GB_readiness/", created_files),
    panel_or_figure = "NAR G&B readiness package",
    change_type = "Priority6_NAR_GB_readiness",
    old_text = "",
    new_text = "Created reproducibility, data/code availability, provenance, manifest, session-info, release-plan, review-inventory, and gap-audit materials.",
    plot_visible_or_caption = "documentation/manifest",
    regenerated_or_copied = "generated",
    notes = "Documentation-only package; no IMRS scoring, validation, robustness, enrichment, figure, or supplementary table values were changed."
  )
  readr::write_tsv(bind_rows(existing[, cols], rows), path, na = "")
}

update_table_inventory <- function() {
  path <- paths$table_inventory
  cols <- c("table_file", "purpose_inferred", "edited_yes_no", "key_changes", "row_count", "column_count")
  existing <- if (file.exists(path)) read_tsv_chr(path) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% filter(!str_starts(.data$table_file, "NAR_GB_readiness/"))
  table_files <- list.files(out_dir, pattern = "\\.tsv$", full.names = TRUE)
  rows <- map_dfr(sort(table_files), function(p) {
    tbl <- read_tsv_chr(p)
    tibble(
      table_file = paste0("NAR_GB_readiness/", basename(p)),
      purpose_inferred = case_when(
        str_detect(basename(p), "accession") ~ "NAR G&B public dataset accession table",
        str_detect(basename(p), "script_to_output") ~ "NAR G&B script-to-output manifest",
        str_detect(basename(p), "input_data") ~ "NAR G&B input data manifest",
        str_detect(basename(p), "derived_output") ~ "NAR G&B derived output manifest",
        str_detect(basename(p), "review_package") ~ "NAR G&B review package inventory",
        str_detect(basename(p), "gap_audit") ~ "NAR G&B submission gap audit",
        str_detect(basename(p), "checklist") ~ "NAR G&B reproducibility checklist",
        TRUE ~ "NAR G&B readiness table"
      ),
      edited_yes_no = "yes",
      key_changes = "Created for Priority 6 NAR G&B readiness package",
      row_count = as.character(nrow(tbl)),
      column_count = as.character(ncol(tbl))
    )
  })
  readr::write_tsv(bind_rows(existing[, cols], rows), path, na = "")
}

update_file_inventory <- function() {
  path <- paths$file_inventory
  cols <- c("file", "extension", "size_bytes", "modified_time", "copied_from_v5_yes_no", "generated_in_v6_yes_no")
  existing <- if (file.exists(path)) read_tsv_chr(path) else tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>% filter(!str_starts(.data$file, "NAR_GB_readiness/"))
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  rows <- map_dfr(sort(files), function(p) {
    info <- file.info(p)
    tibble(
      file = paste0("NAR_GB_readiness/", basename(p)),
      extension = tools::file_ext(p),
      size_bytes = as.character(as.numeric(info$size)),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      copied_from_v5_yes_no = "no",
      generated_in_v6_yes_no = "yes"
    )
  })
  readr::write_tsv(bind_rows(existing[, cols], rows), path, na = "")
}

update_change_log()
update_table_inventory()

required_files <- c(
  "IMRS_reproducibility_README_v6.md",
  "IMRS_code_availability_statement_candidate.txt",
  "IMRS_data_availability_statement_candidate.txt",
  "IMRS_reproducibility_checklist_v6.tsv",
  "IMRS_public_dataset_accession_table_v6.tsv",
  "IMRS_script_to_output_manifest_v6.tsv",
  "IMRS_input_data_manifest_v6.tsv",
  "IMRS_derived_output_manifest_v6.tsv",
  "IMRS_software_session_info_v6.txt",
  "IMRS_repository_release_plan_v6.md",
  "IMRS_review_package_inventory_v6.tsv",
  "IMRS_NAR_GB_submission_gap_audit_v6.tsv",
  "Priority6_NAR_GB_readiness_log.txt"
)
required_paths <- file.path(out_dir, required_files)
validation_lines <- c(
  paste0("NAR_GB_readiness folder exists: ", dir.exists(out_dir)),
  paste0("Required files created and nonzero before final log write, except log: ", all(file.exists(required_paths[required_files != "Priority6_NAR_GB_readiness_log.txt"]) & file.info(required_paths[required_files != "Priority6_NAR_GB_readiness_log.txt"])$size > 0)),
  paste0("Public dataset accession rows: ", nrow(accession_table)),
  paste0("Script-to-output manifest rows: ", nrow(script_manifest)),
  paste0("Input data manifest rows: ", nrow(input_manifest)),
  paste0("Derived output manifest rows: ", nrow(derived_manifest)),
  paste0("Review package inventory rows: ", nrow(review_inventory)),
  paste0("NAR G&B gap audit rows: ", nrow(gap_audit)),
  "R sessionInfo captured in IMRS_software_session_info_v6.txt.",
  "No IMRS scoring, validation, robustness, enrichment, figure, or supplementary table values were modified.",
  "No new mechanistic, clinical prediction, or delivery-platform safety-ranking claims were introduced."
)
log_msg("Validation summary:")
for (line in validation_lines) log_msg("  ", line)
log_msg("Manual actions remaining: confirm raw-data redistribution/retrieval policy and add an archival DOI only if one is issued in a future deposit.")
log_msg("Updated v6_change_log.tsv, v6_file_inventory.tsv, and v6_table_inventory.tsv.")
log_msg("Priority 6 NAR G&B readiness packaging complete.")
write_log()
update_file_inventory()

cat("Priority 6 NAR G&B readiness package complete\n")
cat("Output folder: ", out_dir, "\n", sep = "")
cat("Public dataset accession rows: ", nrow(accession_table), "\n", sep = "")
cat("Script-to-output manifest rows: ", nrow(script_manifest), "\n", sep = "")
cat("Input data manifest rows: ", nrow(input_manifest), "\n", sep = "")
cat("Derived output manifest rows: ", nrow(derived_manifest), "\n", sep = "")
cat("Warnings: ", ifelse(length(warnings_found) == 0, "none", paste(warnings_found, collapse = "; ")), "\n", sep = "")
