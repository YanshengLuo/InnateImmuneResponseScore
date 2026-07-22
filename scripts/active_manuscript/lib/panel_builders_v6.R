#!/usr/bin/env Rscript

# Master IMRS reorganized figure-generation script.
# Regenerates current manuscript/PPT figure panels from existing result tables only.
# Outputs are written to revised_plots and never overwrite manuscript source
# figure directories.

options(stringsAsFactors = FALSE)

# -------------------------------------------------------------------------
# A. Setup
# -------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble", "purrr",
  "ggplot2", "scales", "grid"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(grid)
})

project_root <- normalizePath(Sys.getenv("IMRS_REPOSITORY_ROOT", "."),
                              winslash = "/", mustWork = FALSE)
figure_input_dir <- normalizePath(
  Sys.getenv("IMRS_FIGURE_INPUT_DIR", file.path(project_root, "data", "derived", "figure_inputs")),
  winslash = "/", mustWork = FALSE
)
if (!dir.exists(figure_input_dir)) {
  stop("Released figure-input directory does not exist: ", figure_input_dir, call. = FALSE)
}
output_root <- normalizePath(
  Sys.getenv("IMRS_MANUSCRIPT_OUTPUT_DIR", file.path(project_root, "results_release_templates", "figures")),
  winslash = "/", mustWork = FALSE
)
final_root <- figure_input_dir
v2_root <- figure_input_dir
scripts_root <- file.path(project_root, "scripts", "active_manuscript", "lib")
existing_figures_root <- output_root
audit_results_dir <- figure_input_dir
extra_results_dir <- figure_input_dir
temporal_root <- file.path(output_root, "_panel_builder_work")
dir.create(temporal_root, recursive = TRUE, showWarnings = FALSE)
five_anchor_loao_root <- file.path(temporal_root, "five_anchor_LOAO")
dir.create(five_anchor_loao_root, recursive = TRUE, showWarnings = FALSE)
figure1c_inclusion_check_path <- file.path(output_root, "Figure1C_pseudolog_inclusion_check.tsv")
figure8a_loao_source_audit_path <- file.path(output_root, "Figure8A_LOAO_source_audit.tsv")
figure8b_sample_source_check_path <- file.path(output_root, "Figure8B_anchor_sample_source_check.tsv")
figure4a_anchor_inclusion_check_path <- file.path(output_root, "Figure4A_anchor_inclusion_check.tsv")

figure_folders <- c(
  "Figure1_framework_dataset_provenance",
  "Figure2_anchor_construction_weights",
  "Figure3_primary_validation_overview",
  "Figure4_validation_detail_and_discrimination",
  "Figure5_permutation_null_analysis",
  "Figure6_baseline_signature_benchmarking",
  "Figure7_gene_threshold_robustness",
  "Figure8_anchor_coefficient_robustness",
  "FigureS1_weak_late_context_summary",
  "FigureS2_detailed_validation_forests",
  "FigureS3_context_timecourse_and_dominance_appendix"
)
folder_path <- function(folder) file.path(output_root, folder)
invisible(lapply(file.path(output_root, figure_folders), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

manifest_path <- file.path(output_root, "figure_generation_manifest.tsv")
log_path <- file.path(output_root, "figure_generation_log.txt")
context_tsv_path <- file.path(output_root, "figure_dataset_context_companion.tsv")
context_md_path <- file.path(output_root, "figure_dataset_context_companion.md")
figure2b_support_pattern_path <- file.path(output_root, "Figure2B_core_gene_support_pattern_notes.tsv")

log_lines <- character()
warnings_seen <- character()
skipped_panels <- tibble(panel_id = character(), reason = character())
manifest_rows <- list()
output_registry <- character()

log_msg <- function(..., level = "INFO") {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 level, " ", paste(..., collapse = ""))
  message(line)
  log_lines <<- c(log_lines, line)
}

add_warning <- function(...) {
  msg <- paste(..., collapse = "")
  warnings_seen <<- unique(c(warnings_seen, msg))
  log_msg(msg, level = "WARN")
}

write_log <- function() {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  cat(paste(log_lines, collapse = "\n"), "\n", file = log_path)
}

locked_anchor_ids <- c(
  "GSE39129",
  "GSE167521",
  "GSE264344",
  "GSE279372",
  "GSE279744"
)
# Strict-3 is a sensitivity/ablation set only, not the production model set.
strict3_sensitivity_ids <- c("GSE39129", "GSE167521", "GSE264344")
anchor_palette <- c(
  GSE39129 = "#4E79A7",
  GSE167521 = "#59A14F",
  GSE264344 = "#E15759",
  GSE279372 = "#7B5EA7",
  GSE279744 = "#F28E2B"
)
strict_anchor_ids <- locked_anchor_ids
LOCKED_DATASETS_MOUSE <- locked_anchor_ids
additional_discovery_palette <- setNames(character(), character())

condition_palette <- c(
  CONTROL = "#6B7280",
  DELIVERY = "#B22222"
)
condition_labels <- c(CONTROL = "Control", DELIVERY = "Delivery")

manuscript_plot_groups <- c(
  "Locked anchor",
  "Primary acute validation",
  "Extended validation",
  "Secondary support"
)
validation_plot_groups <- setdiff(manuscript_plot_groups, "Locked anchor")
manuscript_group_order <- manuscript_plot_groups
excluded_unclear_role_values <- c(
  "Excluded/unclear",
  "Excluded",
  "Unclear",
  "excluded/unclear",
  "excluded",
  "unclear",
  "excluded_or_unclear"
)
manuscript_group_palette <- c(
  "Locked anchor" = "#2F6B9A",
  "Primary acute validation" = "#C84B31",
  "Extended validation" = "#7B5EA7",
  "Secondary support" = "#59A14F"
)
dataset_role_summary_labels <- c(
  "Locked anchor" = "Exact anchor",
  "Primary acute validation" = "Primary acute validation",
  "Extended validation" = "Extended validation",
  "Secondary support" = "Secondary support"
)

score_palette <- c(
  "IMRS" = "#111827",
  "ISG baseline" = "#0072B2",
  "Inflammatory baseline" = "#D55E00",
  "Other benchmark" = "#009E73"
)

reviewer_risk_palette <- c(
  low = "#2F6B9A",
  medium = "#F28E2B",
  high = "#E15759"
)

weight_direction_palette <- c(
  "Positive weight" = "#4E79A7",
  "Negative weight" = "#D55E00"
)

preservation_palette <- c(
  "Direction preserved" = "#4E79A7",
  "Any direction change" = "#B22222"
)

threshold_fdr_palette <- c(
  "0.05" = "#4E79A7",
  "0.1" = "#59A14F",
  "0.2" = "#F28E2B"
)

theme_imrs_publication <- function(base_size = 11, legend_position = "bottom") {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "#111827", size = base_size + 2),
      plot.subtitle = element_text(color = "#374151", size = base_size),
      axis.title = element_text(face = "bold", color = "#111827"),
      axis.text = element_text(color = "#111827"),
      legend.position = legend_position,
      legend.title = element_text(face = "bold", color = "#111827"),
      legend.text = element_text(color = "#111827"),
      strip.text = element_text(face = "bold", color = "#111827"),
      plot.margin = margin(8, 10, 8, 10),
      panel.grid.minor = element_blank()
    )
}

required_paths <- list(
  role_table = file.path(figure_input_dir, "manuscript_dataset_role_table.tsv"),
  dataset_classification = file.path(figure_input_dir, "dataset_classification_authoritative.tsv"),
  leave_one_gene_out = file.path(figure_input_dir, "leave_one_gene_out_summary.tsv"),
  gene_dominance = file.path(figure_input_dir, "gene_dominance_summary.tsv"),
  threshold_summary = file.path(figure_input_dir, "threshold_sensitivity_summary.tsv"),
  threshold_detail = file.path(figure_input_dir, "threshold_sensitivity_contrast_deltas.tsv"),
  leave_one_anchor_summary = file.path(figure_input_dir, "five_anchor_leave_one_anchor_out_summary.tsv"),
  leave_one_anchor_detail = file.path(figure_input_dir, "five_anchor_leave_one_anchor_out_contrast_details.tsv"),
  weak_context = file.path(figure_input_dir, "weak_dataset_paper_context_audit.tsv"),
  anchor_support = file.path(figure_input_dir, "support_by_dataset.tsv"),
  anchor_contrast_counts_06A = file.path(figure_input_dir, "contrast_counts_by_dataset.tsv"),
  anchor_contrast_counts_06B = file.path(figure_input_dir, "contrast_counts_by_dataset_6B.tsv"),
  gene_weights = file.path(figure_input_dir, "gene_weights.tsv"),
  gene_symbols = file.path(figure_input_dir, "gene_symbol_mapping.tsv"),
  gene_heterogeneity = file.path(figure_input_dir, "gene_heterogeneity.tsv"),
  gene_power = file.path(figure_input_dir, "gene_power.tsv"),
  step09_eval = file.path(figure_input_dir, "step09_split_eval.tsv"),
  step09_sample = file.path(figure_input_dir, "step09_split_sample_level.tsv"),
  label_permutation_summary = file.path(figure_input_dir, "label_permutation_null_summary.tsv"),
  baseline_contrast_long = file.path(figure_input_dir, "baseline_signature_contrast_long.tsv"),
  baseline_paired = file.path(figure_input_dir, "baseline_signature_paired_contrast_comparison.tsv"),
  baseline_summary = file.path(figure_input_dir, "baseline_signature_summary_by_group.tsv"),
  baseline_scores_sample = file.path(figure_input_dir, "baseline_signature_scores_sample_level.tsv"),
  coefficient_sensitivity_summary = file.path(figure_input_dir, "coefficient_sensitivity_summary.tsv")
)

check_required_paths <- function(paths) {
  missing <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
  if (length(missing) > 0) {
    missing_text <- paste(paste0(missing, ": ", unlist(paths[missing], use.names = FALSE)), collapse = "\n")
    stop("Missing required input table(s):\n", missing_text, call. = FALSE)
  }
}

read_required_tsv <- function(path, required_cols = character(), label = basename(path)) {
  if (!file.exists(path)) stop("Missing required input table: ", path, call. = FALSE)
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, guess_max = 100000)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Input table ", label, " is missing required column(s): ",
         paste(missing_cols, collapse = ", "), "\nPath: ", path, call. = FALSE)
  }
  df
}

check_required_cols <- function(df, cols, label) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(label, " is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

check_not_all_na <- function(df, col, label) {
  if (!col %in% names(df)) stop(label, " is missing required column: ", col, call. = FALSE)
  if (all(is.na(df[[col]]))) stop(label, " has NA-only plotting variable: ", col, call. = FALSE)
  invisible(TRUE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))
safe_chr <- function(x) ifelse(is.na(x) | !nzchar(as.character(x)), "Not specified", as.character(x))
logic_col <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}
strip_ens <- function(x) sub("\\.\\d+$", "", as.character(x))
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) hit[[1]] else NA_character_
}
display_text <- function(x) stringr::str_replace_all(as.character(x), "_", " ")
short_text <- function(x, n = 64) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
collapse_paths <- function(x) paste(unique(as.character(x)), collapse = ";")
ordered_factor <- function(label, order_value, decreasing = FALSE) {
  label <- as.character(label)
  order_value <- safe_num(order_value)
  ord <- order(order_value, decreasing = decreasing, na.last = TRUE)
  factor(label, levels = unique(label[ord]))
}

format_time_label <- function(time_h) {
  time_h <- safe_num(time_h)
  ifelse(is.finite(time_h), paste0(format(time_h, trim = TRUE, scientific = FALSE), " h"), "time not specified")
}

dataset_context_one <- function(dataset_id, tissue = NA_character_, time_h = NA_real_,
                                platform = NA_character_, compact = FALSE) {
  dataset_id <- as.character(dataset_id)
  tissue <- safe_chr(tissue)
  platform <- safe_chr(platform)
  time_h_num <- safe_num(time_h)
  time_label <- format_time_label(time_h_num)
  tissue_label <- display_text(tissue)
  platform_label <- display_text(platform)
  if (dataset_id == "GSE39129") return("GSE39129 | lentiviral vector | liver | 4 h")
  if (dataset_id == "GSE167521") return("GSE167521 | LNP | skin | acute")
  if (dataset_id == "GSE264344") {
    if (is.finite(time_h_num) && time_h_num <= 24) {
      return(paste("GSE264344 | adenoviral vector |", tissue_label, "| 1-24 h"))
    }
    if (is.finite(time_h_num) && time_h_num > 24) {
      return(paste("GSE264344 | adenoviral vector |", tissue_label, "|", time_label))
    }
    return("GSE264344 | adenoviral vector | muscle/dLN/blood | time course")
  }
  if (dataset_id == "GSE279744") return("GSE279744 | LNP formulation | dLN")
  if (compact) return(paste(dataset_id, tissue_label, time_label, sep = " | "))
  paste(dataset_id, platform_label, tissue_label, time_label, sep = " | ")
}

dataset_context_label <- function(dataset_id, tissue = NA_character_, time_h = NA_real_,
                                  platform = NA_character_, compact = FALSE) {
  mapply(dataset_context_one, dataset_id, tissue, time_h, platform,
         MoreArgs = list(compact = compact), USE.NAMES = FALSE)
}

normalize_manuscript_group <- function(code_or_label, fallback_role = NA_character_) {
  x <- as.character(code_or_label)
  role <- as.character(fallback_role)
  out <- dplyr::case_when(
    x %in% manuscript_group_order ~ x,
    x %in% c("anchor", "anchor_discovery", "locked_anchor", "strict_platform_anchor", "strict_anchor", "derivation_anchor",
             "additional_acute_discovery_support", "discovery_support") |
      role %in% c("anchor", "anchor_discovery", "locked_anchor", "strict_platform_anchor", "strict_anchor", "derivation_anchor",
                  "additional_acute_discovery_support", "discovery_support") ~ "Locked anchor",
    x == "primary_acute_validation" | role == "primary_acute_validation" |
      x == "external_acute" | role == "external_acute" |
      x == "primary_external_validation" | role == "primary_external_validation" ~ "Primary acute validation",
    x == "extended_validation" | role == "extended_validation" |
      x == "external_extended" | role == "external_extended" |
      x == "extended_exploratory_transfer" | role == "extended_exploratory_transfer" ~ "Extended validation",
    x == "secondary_support" | role == "secondary_support" |
      x == "secondary_support_not_primary" | role == "secondary_support_not_primary" ~ "Secondary support",
    x == "calibration" ~ "Secondary support",
    x == "external_acute" ~ "Primary acute validation",
    x == "external_extended" ~ "Extended validation",
    TRUE ~ "Excluded/unclear"
  )
  factor(out, levels = manuscript_group_order)
}

is_excluded_unclear_value <- function(x) {
  y <- stringr::str_squish(as.character(x))
  !is.na(y) & nzchar(y) &
    stringr::str_to_lower(y) %in% stringr::str_to_lower(excluded_unclear_role_values)
}

count_excluded_unclear_rows <- function(df, role_cols) {
  role_cols <- intersect(role_cols, names(df))
  if (length(role_cols) == 0) return(0L)
  flags <- lapply(role_cols, function(col) is_excluded_unclear_value(df[[col]]))
  sum(Reduce(`|`, flags), na.rm = TRUE)
}

clean_dataset_role_for_plot <- function(df, role_col = "manuscript_group",
                                        role_levels = manuscript_plot_groups) {
  if (!role_col %in% names(df)) {
    stop("Missing role column for plot cleaning: ", role_col, call. = FALSE)
  }
  role_values <- stringr::str_squish(as.character(df[[role_col]]))
  excluded_rows <- is_excluded_unclear_value(role_values)
  keep_rows <- !is.na(role_values) & nzchar(role_values) & !excluded_rows
  out <- df[keep_rows, , drop = FALSE]
  out[[role_col]] <- factor(
    stringr::str_squish(as.character(out[[role_col]])),
    levels = role_levels
  )
  out <- out[!is.na(out[[role_col]]), , drop = FALSE]
  attr(out, "excluded_unclear_rows_removed") <- sum(excluded_rows, na.rm = TRUE)
  out
}

map_score_label <- function(score_label, score_id = NA_character_) {
  label <- tolower(paste(score_label, score_id))
  out <- dplyr::case_when(
    str_detect(label, "\\bimrs\\b") ~ "IMRS",
    str_detect(label, "isg|interferon|ifit|mx1|oas") ~ "ISG baseline",
    str_detect(label, "inflam|chemokine|ccl|cxcl|il1|il6|tnf") ~ "Inflammatory baseline",
    TRUE ~ "Other benchmark"
  )
  factor(out, levels = names(score_palette))
}

save_imrs_plot <- function(plot, out_dir, stem, width, height, dpi = 400,
                           source_tables = character(),
                           source_code_section_or_function = NA_character_,
                           notes = NA_character_) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  registry_key <- paste(normalizePath(out_dir, winslash = "/", mustWork = FALSE), stem, sep = "::")
  if (registry_key %in% output_registry) {
    stop("Duplicated output stem inside same folder: ", stem, " in ", out_dir, call. = FALSE)
  }
  output_registry <<- c(output_registry, registry_key)

  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  plot_notes <- notes

  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi,
                  limitsize = FALSE, bg = "white")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, limitsize = FALSE, bg = "white")

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(svg_path, plot, width = width, height = height,
                    device = svglite::svglite, limitsize = FALSE, bg = "white")
  } else {
    svg_path <- NA_character_
    plot_notes <- paste(c(plot_notes, "SVG skipped because svglite is not installed."), collapse = " ")
    add_warning("SVG skipped for ", stem, " because svglite is not installed.")
  }

  row <- tibble(
    figure_id = basename(out_dir),
    panel_id = stem,
    output_folder = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    output_png = normalizePath(png_path, winslash = "/", mustWork = FALSE),
    output_pdf = normalizePath(pdf_path, winslash = "/", mustWork = FALSE),
    output_svg = ifelse(is.na(svg_path), NA_character_,
                        normalizePath(svg_path, winslash = "/", mustWork = FALSE)),
    source_tables = collapse_paths(source_tables),
    source_code_section_or_function = source_code_section_or_function,
    width = width,
    height = height,
    dpi = dpi,
    notes = plot_notes
  )
  manifest_rows[[length(manifest_rows) + 1L]] <<- row
  log_msg("Saved ", stem, " to ", out_dir)
  invisible(row)
}

save_imrs_grid <- function(draw_fun, out_dir, stem, width, height, dpi = 400,
                           source_tables = character(),
                           source_code_section_or_function = NA_character_,
                           notes = NA_character_) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  registry_key <- paste(normalizePath(out_dir, winslash = "/", mustWork = FALSE), stem, sep = "::")
  if (registry_key %in% output_registry) {
    stop("Duplicated output stem inside same folder: ", stem, " in ", out_dir, call. = FALSE)
  }
  output_registry <<- c(output_registry, registry_key)

  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  grid_notes <- notes

  grDevices::png(png_path, width = width, height = height, units = "in",
                 res = dpi, bg = "white")
  draw_fun()
  grDevices::dev.off()

  grDevices::cairo_pdf(pdf_path, width = width, height = height, onefile = FALSE, bg = "white")
  draw_fun()
  grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_path, width = width, height = height, bg = "white")
    draw_fun()
    grDevices::dev.off()
  } else {
    svg_path <- NA_character_
    grid_notes <- paste(c(grid_notes, "SVG skipped because svglite is not installed."), collapse = " ")
    add_warning("SVG skipped for grid panel ", stem, " because svglite is not installed.")
  }

  row <- tibble(
    figure_id = basename(out_dir),
    panel_id = stem,
    output_folder = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    output_png = normalizePath(png_path, winslash = "/", mustWork = FALSE),
    output_pdf = normalizePath(pdf_path, winslash = "/", mustWork = FALSE),
    output_svg = ifelse(is.na(svg_path), NA_character_,
                        normalizePath(svg_path, winslash = "/", mustWork = FALSE)),
    source_tables = collapse_paths(source_tables),
    source_code_section_or_function = source_code_section_or_function,
    width = width,
    height = height,
    dpi = dpi,
    notes = grid_notes
  )
  manifest_rows[[length(manifest_rows) + 1L]] <<- row
  log_msg("Saved grid panel ", stem, " to ", out_dir)
  invisible(row)
}

validate_png_outputs <- function(manifest) {
  bad <- manifest %>%
    filter(is.na(output_png) | !file.exists(output_png) | file.info(output_png)$size <= 0)
  if (nrow(bad) > 0) {
    stop("One or more generated PNG files are missing or zero-byte: ",
         paste(bad$panel_id, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_manifest <- function(manifest) {
  dup <- manifest %>% count(output_folder, panel_id) %>% filter(n > 1)
  if (nrow(dup) > 0) {
    stop("Duplicated output stem(s) in manifest: ", paste(dup$panel_id, collapse = ", "), call. = FALSE)
  }
  validate_png_outputs(manifest)
  expected_anchor_palette <- c(
    GSE39129 = "#4E79A7",
    GSE167521 = "#59A14F",
    GSE264344 = "#E15759",
    GSE279372 = "#7B5EA7",
    GSE279744 = "#F28E2B"
  )
  if (!identical(anchor_palette[names(expected_anchor_palette)], expected_anchor_palette) ||
      !identical(names(anchor_palette), names(expected_anchor_palette))) {
    stop("Anchor palette failed consistency check.", call. = FALSE)
  }
  if (!identical(unname(condition_palette["CONTROL"]), "#6B7280") ||
      !identical(unname(condition_palette["DELIVERY"]), "#B22222")) {
    stop("CONTROL/DELIVERY palette failed consistency check.", call. = FALSE)
  }
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# B. Load shared input tables
# -------------------------------------------------------------------------

log_msg("Starting IMRS figure generation.")
log_msg("Project root: ", project_root)
log_msg("Existing figures root will not be modified: ", existing_figures_root)
log_msg("Regenerated output root: ", output_root)
log_msg("Standalone panel mode: no combined A-I figure images are generated.")

check_required_paths(required_paths)

role_tbl <- read_required_tsv(
  required_paths$role_table,
  c("dataset_id", "gse_id", "split_id", "manuscript_role", "final_display_group_v2",
    "tissue", "time_h",
    "delivery_platform_clean", "n_controls", "n_delivery", "delta_mean_imrs_z",
    "auc_imrs_z", "pass"),
  "manuscript_dataset_role_table.tsv"
) %>%
  mutate(
    pass = logic_col(pass),
    time_h = safe_num(time_h),
    delta_mean_imrs_z = safe_num(delta_mean_imrs_z),
    auc_imrs_z = safe_num(auc_imrs_z),
    n_controls = safe_num(n_controls),
    n_delivery = safe_num(n_delivery),
    locked_anchor_dataset = dataset_id %in% LOCKED_DATASETS_MOUSE,
    row_is_anchor_phase = case_when(
      !is.na(time_h) & time_h <= 24 ~ TRUE,
      !is.na(time_h) & time_h > 24 ~ FALSE,
      TRUE ~ locked_anchor_dataset
    ),
    production_anchor_row = locked_anchor_dataset & row_is_anchor_phase,
    final_display_group_v2 = case_when(
      production_anchor_row ~ "Locked anchor",
      dataset_id == "GSE264344" & !is.na(time_h) & time_h > 24 ~ "Extended validation",
      str_detect(dataset_id, "GSE262515") ~ "Secondary support",
      dataset_id %in% c("GSE119119", "GSE139529", "GSE279743") ~ "Primary acute validation",
      dataset_id %in% c("GSE166655", "GSE178313", "GSE314070") ~ "Extended validation",
      TRUE ~ as.character(final_display_group_v2)
    ),
    role_display = final_display_group_v2,
    manuscript_role = case_when(
      production_anchor_row ~ "anchor",
      dataset_id == "GSE264344" & !is.na(time_h) & time_h > 24 ~ "extended_validation",
      str_detect(dataset_id, "GSE262515") ~ "secondary_support",
      dataset_id %in% c("GSE119119", "GSE139529", "GSE279743") ~ "primary_acute_validation",
      dataset_id %in% c("GSE166655", "GSE178313", "GSE314070") ~ "extended_validation",
      TRUE ~ as.character(manuscript_role)
    ),
    manuscript_claim_group = case_when(
      production_anchor_row ~ "anchor_discovery",
      manuscript_role == "primary_acute_validation" ~ "primary_external_validation",
      manuscript_role == "extended_validation" ~ "extended_exploratory_transfer",
      manuscript_role == "secondary_support" ~ "secondary_support_not_primary",
      TRUE ~ as.character(manuscript_claim_group)
    ),
    manuscript_group = factor(final_display_group_v2, levels = manuscript_group_order),
    context_label = dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean)
  )
role_pass <- role_tbl %>% filter(pass %in% TRUE)
role_pass_excluded_unclear_rows_removed <- count_excluded_unclear_rows(
  role_pass,
  c("final_display_group_v2", "manuscript_group", "manuscript_role",
    "manuscript_claim_group", "manuscript_interpretation_group", "role_display")
)
role_pass_for_plot <- clean_dataset_role_for_plot(role_pass, "manuscript_group")
check_not_all_na(role_pass, "delta_mean_imrs_z", "role_pass")

role_lookup <- role_pass_for_plot %>%
  select(gse_id, dataset_id, split_id, manuscript_group, context_label, delivery_platform_clean) %>%
  distinct()

attach_role_group <- function(tbl) {
  out <- tbl
  if (all(c("gse_id", "split_id") %in% names(out))) {
    out <- out %>%
      left_join(
        role_lookup %>% select(gse_id, split_id, manuscript_group_v2 = manuscript_group),
        by = c("gse_id", "split_id")
      )
  } else if (all(c("dataset_id", "split_id") %in% names(out))) {
    out <- out %>%
      left_join(
        role_lookup %>% select(dataset_id, split_id, manuscript_group_v2 = manuscript_group),
        by = c("dataset_id", "split_id")
      )
  } else {
    out$manuscript_group_v2 <- NA_character_
  }
  fallback <- if ("final_display_group_v2" %in% names(out)) {
    out$final_display_group_v2
  } else if ("manuscript_interpretation_group" %in% names(out)) {
    out$manuscript_interpretation_group
  } else if ("manuscript_role" %in% names(out)) {
    out$manuscript_role
  } else {
    NA_character_
  }
  role <- if ("manuscript_role" %in% names(out)) out$manuscript_role else NA_character_
  out <- out %>%
    mutate(
      manuscript_group = factor(
        ifelse(!is.na(manuscript_group_v2),
               as.character(manuscript_group_v2),
               as.character(normalize_manuscript_group(fallback, role))),
        levels = manuscript_group_order
      )
    ) %>%
    select(-manuscript_group_v2)
  clean_dataset_role_for_plot(out, "manuscript_group")
}

dataset_classification_tbl <- read_required_tsv(
  required_paths$dataset_classification,
  c("dataset_id", "used_in_step06A_core_construction", "used_in_step06B_heterogeneity",
    "strict_platform_anchor", "discovery_support_role", "final_display_group"),
  "dataset_classification_authoritative.tsv"
) %>%
  mutate(
    used_in_step06A_core_construction = logic_col(used_in_step06A_core_construction),
    used_in_step06B_heterogeneity = logic_col(used_in_step06B_heterogeneity),
    strict_platform_anchor = logic_col(strict_platform_anchor)
  )

support_tbl <- read_required_tsv(required_paths$anchor_support,
                                 c("dataset_id", "gene_id", "dataset_support_flag"),
                                 "support_by_dataset.tsv")
contrast_counts_06A_tbl <- read_required_tsv(required_paths$anchor_contrast_counts_06A,
                                             c("dataset_id", "n_contrasts"),
                                             "contrast_counts_by_dataset.tsv")
contrast_counts_06B_tbl <- read_required_tsv(required_paths$anchor_contrast_counts_06B,
                                             c("dataset_id", "n_contrasts"),
                                             "contrast_counts_by_dataset_6B.tsv")
weights_tbl <- read_required_tsv(required_paths$gene_weights, character(), "gene_weights.tsv")
symbols_tbl <- read_required_tsv(required_paths$gene_symbols, c("ensembl_gene_id", "mgi_symbol"),
                                 "gene_symbol_mapping.tsv")
heterogeneity_tbl <- read_required_tsv(required_paths$gene_heterogeneity,
                                       c("gene_id", "heterogeneity_flag"),
                                       "gene_heterogeneity.tsv")
power_tbl <- read_required_tsv(required_paths$gene_power,
                               c("gene_id", "low_power_flag"),
                               "gene_power.tsv")

loo_tbl <- read_required_tsv(required_paths$leave_one_gene_out,
                             c("removed_gene_id", "removed_gene_symbol", "gse_id", "split_id",
                               "original_delta_mean_imrs_z", "leave_one_gene_out_delta_mean_imrs_z",
                               "absolute_change_delta", "direction_preserved"),
                             "leave_one_gene_out_summary.tsv")
dominance_tbl <- read_required_tsv(required_paths$gene_dominance,
                                   c("gse_id", "split_id", "median_max_contribution_fraction"),
                                   "gene_dominance_summary.tsv")
threshold_summary_tbl <- read_required_tsv(required_paths$threshold_summary,
                                           c("min_abs_log2FC", "min_up_anchor_count", "fdr_support",
                                             "n_genes_selected", "external_mean_delta_imrs_z"),
                                           "threshold_sensitivity_summary.tsv")
threshold_detail_tbl <- read_required_tsv(required_paths$threshold_detail,
                                          c("min_abs_log2FC", "min_up_anchor_count", "fdr_support",
                                            "delta_mean_imrs_z", "sensitivity_scope"),
                                          "threshold_sensitivity_contrast_deltas.tsv")
loao_summary_tbl <- read_required_tsv(required_paths$leave_one_anchor_summary,
                                      c("holdout_anchor"),
                                      "leave_one_anchor_out_summary.tsv")
loao_detail_tbl <- read_required_tsv(required_paths$leave_one_anchor_detail,
                                     c("holdout_anchor", "delta_mean_imrs_z"),
                                     "leave_one_anchor_out_contrast_details.tsv")
weak_tbl <- read_required_tsv(required_paths$weak_context,
                              c("dataset_id", "gse_id", "split_id", "tissue", "time_h",
                                "treatment_group", "delivery_platform", "original_IMRS_delta",
                                "explanation_category", "reviewer_risk_level"),
                              "weak_dataset_paper_context_audit.tsv")
step09_eval_tbl <- read_required_tsv(required_paths$step09_eval,
                                     c("gse_id", "dataset_id", "split_id", "tissue", "time_h",
                                       "delta_mean_imrs_z", "auc_imrs_z", "pass"),
                                     "step09_split_eval.tsv")
step09_sample_tbl <- read_required_tsv(required_paths$step09_sample,
                                       c("gse_id", "dataset_id", "split_id", "sample_id",
                                         "condition_simple", "imrs_z"),
                                       "step09_split_sample_level.tsv")

perm_summary_tbl <- read_required_tsv(required_paths$label_permutation_summary,
                                      c("gse_id", "split_id", "manuscript_interpretation_group",
                                        "observed_delta_mean_imrs_z", "null_q025_delta",
                                        "null_q975_delta", "empirical_p_greater"),
                                      "label_permutation_null_summary.tsv")
baseline_long_tbl <- read_required_tsv(required_paths$baseline_contrast_long,
                                       c("score_id", "score_label", "manuscript_interpretation_group",
                                         "delta_score", "pass"),
                                       "baseline_signature_contrast_long.tsv")
baseline_paired_tbl <- read_required_tsv(required_paths$baseline_paired,
                                         c("signature_id", "signature_label", "baseline_delta_score",
                                           "imrs_delta_score"),
                                         "baseline_signature_paired_contrast_comparison.tsv")
baseline_summary_tbl <- read_required_tsv(required_paths$baseline_summary,
                                          c("score_id", "score_label", "manuscript_interpretation_group",
                                            "proportion_positive_delta"),
                                          "baseline_signature_summary_by_group.tsv")
baseline_scores_tbl <- read_required_tsv(required_paths$baseline_scores_sample,
                                         c("dataset_id", "sample_id", "signature_id", "signature_label",
                                           "signature_z"),
                                         "baseline_signature_scores_sample_level.tsv")
coefficient_summary_tbl <- read_required_tsv(required_paths$coefficient_sensitivity_summary,
                                             c("metric", "value", "interpretation", "source_file"),
                                             "coefficient_sensitivity_summary.tsv")

preferred_discovery_order <- c(strict_anchor_ids, names(additional_discovery_palette))
actual_discovery_ids <- dataset_classification_tbl %>%
  filter(used_in_step06A_core_construction %in% TRUE) %>%
  pull(dataset_id) %>%
  as.character() %>%
  unique()
discovery_dataset_ids <- c(
  preferred_discovery_order[preferred_discovery_order %in% actual_discovery_ids],
  setdiff(actual_discovery_ids, preferred_discovery_order)
)
actual_discovery_k <- length(discovery_dataset_ids)
discovery_support_threshold <- ceiling((2 / 3) * actual_discovery_k)
if (actual_discovery_k == 0) {
  stop("No Step 06A anchor-phase discovery datasets were detected in contrast_counts_by_dataset.tsv.", call. = FALSE)
}
if (!setequal(discovery_dataset_ids, LOCKED_DATASETS_MOUSE)) {
  stop("Anchor-facing plots require exactly the five locked anchors: ",
       paste(LOCKED_DATASETS_MOUSE, collapse = ", "),
       ". Detected Step 06A discovery set: ",
       paste(discovery_dataset_ids, collapse = ", "), call. = FALSE)
}
if ("GSE262515" %in% discovery_dataset_ids) {
  add_warning("GSE262515 appears in Step 06A discovery counts; manuscript wording requires manual review.")
}

make_discovery_palette <- function(ids) {
  pal <- c(anchor_palette, additional_discovery_palette)
  missing_ids <- setdiff(ids, names(pal))
  if (length(missing_ids) > 0) {
    pal <- c(pal, setNames(rep("#9CA3AF", length(missing_ids)), missing_ids))
  }
  pal[ids]
}
discovery_palette <- make_discovery_palette(discovery_dataset_ids)
discovery_support_label <- paste0(
  "support in at least ", discovery_support_threshold, " of ",
  actual_discovery_k, " discovery datasets"
)
discovery_set_phrase <- if (identical(sort(discovery_dataset_ids), sort(LOCKED_DATASETS_MOUSE))) {
  "five locked acute mouse anchor datasets"
} else if (identical(sort(discovery_dataset_ids), sort(strict_anchor_ids))) {
  "strict-three-anchor sensitivity discovery set"
} else if ("GSE262515" %in% discovery_dataset_ids) {
  "expanded locked early mouse discovery set"
} else {
  paste0(actual_discovery_k, "-dataset locked acute mouse discovery set")
}
log_msg("Detected Step 06A discovery set: ", paste(discovery_dataset_ids, collapse = ", "),
        " (K=", actual_discovery_k, "; threshold=", discovery_support_threshold, ").")

write_dataset_context_companion <- function() {
  context_tbl <- role_pass_for_plot %>%
    transmute(
      dataset_id,
      gse_id,
      tissue,
      time_h,
      delivery_platform = delivery_platform_clean,
      manuscript_group = as.character(manuscript_group),
      plot_context_label = context_label,
      biological_context_note = case_when(
        dataset_id == "GSE39129" ~ "Lentiviral-vector liver dataset used as a locked anchor.",
        dataset_id == "GSE167521" ~ "LNP skin acute-response dataset used as a locked anchor.",
        dataset_id == "GSE264344" & time_h <= 24 ~ "Adenoviral-vector early time-course locked-anchor context across muscle, draining lymph node, and blood.",
        dataset_id == "GSE264344" & time_h > 24 ~ "Adenoviral-vector late time-course context; weaker 72 h signals are interpreted as waning kinetics.",
        dataset_id == "GSE279372" ~ "Acute mouse discovery context included in the five-dataset locked anchor set.",
        dataset_id == "GSE279744" ~ "LNP formulation context in draining lymph node included in the five-dataset locked anchor set.",
        TRUE ~ "See manuscript role table for full tissue, time, and treatment/control context."
      )
    ) %>%
    distinct()
  readr::write_tsv(context_tbl, context_tsv_path, na = "NA")
  md_lines <- c(
    "# IMRS Figure Dataset Context Companion",
    "",
    "Use this companion when a plot keeps compact GSE-heavy labels for readability.",
    "",
    "| Dataset | Plot context label | Biological context note |",
    "| --- | --- | --- |"
  )
  md_rows <- context_tbl %>%
    distinct(dataset_id, plot_context_label, biological_context_note) %>%
    arrange(dataset_id, plot_context_label) %>%
    mutate(row = paste0("| ", dataset_id, " | ", plot_context_label, " | ",
                        biological_context_note, " |")) %>%
    pull(row)
  cat(paste(c(md_lines, md_rows), collapse = "\n"), "\n", file = context_md_path)
  log_msg("Wrote dataset context companion files.")
}

# -------------------------------------------------------------------------
# C. Plot functions
# -------------------------------------------------------------------------

make_Figure1A <- function() {
  draw_fun <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text("IMRS computation workflow", x = 0.5, y = 0.965,
              gp = gpar(fontsize = 16, fontface = "bold", col = "#111827"))
    grid.text("Locked acute discovery construction, frozen-weight scoring, and split-level transfer evaluation",
              x = 0.5, y = 0.932, gp = gpar(fontsize = 9.5, col = "#374151"))
    boxes <- tibble(
      label = c(
        "Raw count matrices +\nverified metadata",
        "Delivery-versus-control\ncontrast definitions",
        "Discovery-set differential-\nexpression evidence",
        "Reproducible acute\nresponse genes",
        "Heterogeneity and\nlow-power gene filters",
        "Frozen discovery-derived\ngene weights",
        "Target dataset\nnormalization",
        "Control-referenced\ngene z-scores",
        "Weighted sample-level\nIMRS score",
        "Control-standardized\nIMRS z-score",
        "Mean delivery-minus-control\nIMRS z-score",
        "Manuscript dataset groups\nand biological context"
      ),
      group = c("input", "input", rep("anchor", 4), rep("score", 4), "eval", "interpret")
    )
    coords <- tibble(
      x = c(rep(0.25, 6), rep(0.75, 6)),
      y = c(seq(0.82, 0.22, length.out = 6), seq(0.82, 0.22, length.out = 6))
    )
    fills <- c(input = "#EAF2F8", anchor = "#E8F5E9",
               score = "#FFF4E6", eval = "#F3E8FF", interpret = "#F8EAEF")
    box_w <- 0.34
    box_h <- 0.078
    for (i in seq_len(nrow(boxes))) {
      grid.roundrect(coords$x[i], coords$y[i], width = box_w, height = box_h,
                     r = unit(0.01, "npc"),
                     gp = gpar(fill = fills[[boxes$group[i]]], col = "#334155", lwd = 1.1))
      grid.text(boxes$label[i], coords$x[i], coords$y[i],
                gp = gpar(fontsize = 8.4, col = "#111827", lineheight = 0.9))
    }
    draw_arrow <- function(x0, y0, x1, y1) {
      grid.lines(c(x0, x1), c(y0, y1),
                 arrow = arrow(length = unit(0.018, "npc"), type = "closed"),
                 gp = gpar(col = "#475569", lwd = 1.05))
    }
    for (i in 1:5) draw_arrow(coords$x[i], coords$y[i] - box_h / 2,
                              coords$x[i + 1], coords$y[i + 1] + box_h / 2)
    draw_arrow(coords$x[6] + box_w / 2, coords$y[6], coords$x[7] - box_w / 2, coords$y[7])
    for (i in 7:11) draw_arrow(coords$x[i], coords$y[i] - box_h / 2,
                               coords$x[i + 1], coords$y[i + 1] + box_h / 2)
    grid.text("Frozen weights are not refit during validation or transfer.",
              x = 0.5, y = 0.105,
              gp = gpar(fontsize = 9.5, fontface = "italic", col = "#334155"))
  }
  save_imrs_grid(draw_fun, folder_path("Figure1_framework_dataset_provenance"),
                 "Figure1A_imrs_model_computation_workflow_v11", 6, 4.25, dpi = 400,
                 source_tables = character(),
                 source_code_section_or_function = "make_Figure1A")
}

make_Figure1B <- function() {
  draw_fun <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    grid.text("IMRS framework and dataset provenance", x = 0.5, y = 0.95,
              gp = gpar(fontsize = 16, fontface = "bold", col = "#111827"))
    grid.text("Verified split definitions to manuscript-level validation groups",
              x = 0.5, y = 0.915, gp = gpar(fontsize = 9.5, col = "#374151"))
    steps <- tibble(
      label = c(
        "Verified metadata and\nraw count matrices",
        "Delivery-versus-control\nsplit definitions",
        "DESeq2-derived\nacute discovery evidence",
        "Frozen gene weights",
        "Control-standardized\nsample scoring",
        "Split-level transfer\nevaluation",
        "Manuscript curation:\ndataset roles and context",
        "Publication-ready\nfigure panels"
      ),
      group = c("input", "input", "model", "model", "score", "score", "curation", "output")
    )
    x_pos <- rep(c(0.28, 0.72), each = 4)
    y_pos <- rep(seq(0.80, 0.32, length.out = 4), times = 2)
    fills <- c(input = "#EAF2F8", model = "#E8F5E9", score = "#FFF4E6",
               curation = "#F3E8FF", output = "#F8EAEF")
    box_w <- 0.34
    box_h <- 0.092
    for (i in seq_len(nrow(steps))) {
      grid.roundrect(x_pos[i], y_pos[i], width = box_w, height = box_h,
                     r = unit(0.01, "npc"),
                     gp = gpar(fill = fills[[steps$group[i]]], col = "#334155", lwd = 1.1))
      grid.text(steps$label[i], x_pos[i], y_pos[i],
                gp = gpar(fontsize = 8.8, col = "#111827", lineheight = 0.9))
    }
    draw_arrow <- function(x0, y0, x1, y1) {
      grid.lines(c(x0, x1), c(y0, y1),
                 arrow = arrow(length = unit(0.018, "npc"), type = "closed"),
                 gp = gpar(col = "#475569", lwd = 1.05))
    }
    for (i in 1:3) draw_arrow(x_pos[i], y_pos[i] - box_h / 2, x_pos[i + 1], y_pos[i + 1] + box_h / 2)
    draw_arrow(x_pos[4] + box_w / 2, y_pos[4], x_pos[5] - box_w / 2, y_pos[5])
    for (i in 5:7) draw_arrow(x_pos[i], y_pos[i] - box_h / 2, x_pos[i + 1], y_pos[i + 1] + box_h / 2)
    grid.text("Dataset group labels are reader-facing summaries used consistently across figures.",
              x = 0.5, y = 0.145,
              gp = gpar(fontsize = 9.2, fontface = "italic", col = "#334155"))
  }
  save_imrs_grid(draw_fun, folder_path("Figure1_framework_dataset_provenance"),
                 "Figure1B_workflow_overview", 6, 3.95, dpi = 400,
                 source_tables = character(),
                 source_code_section_or_function = "make_Figure1B")
}

make_Figure1C <- function() {
  figure1c_input <- role_pass_for_plot %>%
    filter(is.finite(delta_mean_imrs_z))
  input_contexts <- figure1c_input %>%
    transmute(
      context_key = paste(dataset_id, tissue, time_h, sep = " | ")
    ) %>%
    distinct()

  plot_tbl <- figure1c_input %>%
    group_by(dataset_id, tissue, time_h, delivery_platform_clean, manuscript_group) %>%
    summarise(
      mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
      n_contrasts = n(),
      .groups = "drop"
    ) %>%
    mutate(
      label = dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE),
      context_key = paste(dataset_id, tissue, time_h, sep = " | "),
      label = ordered_factor(label, mean_delta)
    )
  plotted_by_dataset <- plot_tbl %>%
    group_by(dataset_id) %>%
    summarise(
      appears_in_plot = TRUE,
      plotted_group = paste(sort(unique(as.character(manuscript_group))), collapse = "; "),
      plotted_label = paste(sort(unique(as.character(label))), collapse = "; "),
      .groups = "drop"
    )
  inclusion_check <- role_pass_for_plot %>%
    group_by(dataset_id) %>%
    summarise(
      n_role_table_rows = n(),
      n_rows_with_delta = sum(is.finite(delta_mean_imrs_z), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(plotted_by_dataset, by = "dataset_id") %>%
    mutate(
      appears_in_plot = replace_na(appears_in_plot, FALSE),
      plotted_group = ifelse(is.na(plotted_group) | !nzchar(plotted_group), "none", plotted_group),
      plotted_label = ifelse(is.na(plotted_label) | !nzchar(plotted_label), "none", plotted_label),
      reason_if_missing = case_when(
        appears_in_plot ~ "plotted",
        n_rows_with_delta == 0 ~ "no usable delta_mean_imrs_z in role table",
        TRUE ~ "missing despite usable delta_mean_imrs_z"
      )
    ) %>%
    arrange(dataset_id)
  readr::write_tsv(inclusion_check, figure1c_inclusion_check_path, na = "NA")
  missing_with_delta <- inclusion_check %>%
    filter(n_rows_with_delta > 0, !appears_in_plot)
  if (nrow(missing_with_delta) > 0) {
    stop("Figure1C omitted dataset(s) with usable delta_mean_imrs_z: ",
         paste(missing_with_delta$dataset_id, collapse = ", "),
         ". See ", figure1c_inclusion_check_path, call. = FALSE)
  }
  log_msg("Figure1C inclusion check: ", n_distinct(figure1c_input$dataset_id),
          " datasets with usable delta_mean_imrs_z plotted as ", nrow(plot_tbl),
          " dataset/tissue/time/platform contexts.")
  check_not_all_na(plot_tbl, "mean_delta", "Figure1C plot table")
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_point(alpha = 0.92) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 2),
                       breaks = c(-2, -1, 0, 1, 2, 5, 10, 25, 50)) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    scale_size_continuous(range = c(2, 6)) +
    labs(
      title = "IMRS response varies across biological contexts",
      subtitle = "Pseudo-log x-axis preserves resolution near zero while keeping large positive responses visible.",
      x = "Mean delivery-minus-control IMRS z-score (pseudo-log scale)",
      y = "Dataset, biological context, and time window",
      color = "Manuscript analysis group",
      size = "Contrasts"
    ) +
    theme_imrs_publication(base_size = 8)
  save_imrs_plot(p, folder_path("Figure1_framework_dataset_provenance"),
                 "Figure1C_dataset_tissue_pseudolog", 7.8, 6.2, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure1C")
}

make_Figure1D <- function() {
  plot_tbl <- role_pass_for_plot %>%
    group_by(dataset_id, tissue, time_h, delivery_platform_clean, manuscript_group) %>%
    summarise(mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              n_contrasts = n(), .groups = "drop") %>%
    mutate(
      label = paste0(dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE),
                     " (n=", n_contrasts, ")"),
      label = ordered_factor(label, mean_delta)
    )
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_point(size = 2.5, alpha = 0.92) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Manuscript analysis groups preserve response structure",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "Dataset context",
      color = "Manuscript analysis group"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure1_framework_dataset_provenance"),
                 "Figure1D_manuscript_dataset_role_forest", 8.5, 6, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure1D",
                 notes = "Uses the shared manuscript-group palette, not a separate role palette.")
}

make_dataset_role_summary_plot <- function() {
  plot_tbl <- role_pass_for_plot %>%
    group_by(manuscript_group) %>%
    summarise(
      n_contrasts = n(),
      n_datasets = n_distinct(dataset_id),
      .groups = "drop"
    ) %>%
    mutate(manuscript_group = factor(as.character(manuscript_group), levels = manuscript_plot_groups))
  max_y <- max(plot_tbl$n_contrasts, na.rm = TRUE)
  if (!is.finite(max_y) || max_y <= 0) max_y <- 1

  p <- ggplot(plot_tbl, aes(x = manuscript_group, y = n_contrasts, fill = manuscript_group)) +
    geom_col(width = 0.66, alpha = 0.86) +
    geom_text(aes(label = paste0(n_contrasts, " contrasts\n",
                                 n_datasets, " datasets")),
              vjust = -0.15, size = 3.1, lineheight = 0.9) +
    scale_x_discrete(labels = dataset_role_summary_labels, drop = FALSE) +
    scale_fill_manual(values = manuscript_group_palette, breaks = manuscript_plot_groups,
                      labels = dataset_role_summary_labels, drop = FALSE) +
    labs(
      title = "Dataset Role Summary",
      x = NULL,
      y = "Passing split contrasts"
    ) +
    coord_cartesian(ylim = c(0, max_y * 1.22)) +
    theme_imrs_publication(base_size = 9, legend_position = "none") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  save_imrs_plot(p, folder_path("Figure1_framework_dataset_provenance"),
                 "Figure1D_dataset_role_summary", 5.8, 4.2, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_dataset_role_summary_plot",
                 notes = "Regenerated in v2 from the manuscript role table after applying the role-cleaning helper.")
}

make_Figure2A <- function() {
  support <- support_tbl %>%
    mutate(dataset_support_flag = logic_col(dataset_support_flag),
           dataset_id = factor(dataset_id, levels = LOCKED_DATASETS_MOUSE)) %>%
    filter(!is.na(dataset_id))
  plot_tbl <- support %>%
    group_by(dataset_id) %>%
    summarise(supporting_genes = sum(dataset_support_flag, na.rm = TRUE), .groups = "drop") %>%
    right_join(tibble(dataset_id = factor(LOCKED_DATASETS_MOUSE, levels = LOCKED_DATASETS_MOUSE)),
               by = "dataset_id") %>%
    mutate(supporting_genes = replace_na(supporting_genes, 0L))
  p <- ggplot(plot_tbl, aes(x = dataset_id, y = supporting_genes, fill = dataset_id)) +
    geom_col(width = 0.62, show.legend = FALSE) +
    scale_fill_manual(values = anchor_palette[LOCKED_DATASETS_MOUSE], drop = FALSE) +
    labs(
      title = "Anchor-level differential-expression support before reproducibility filtering",
      subtitle = str_wrap(paste0("Shown are the five locked anchors: ",
                                 paste(LOCKED_DATASETS_MOUSE, collapse = ", "),
                                 ". High single-anchor support does not mean that dataset dominates final frozen weights."), width = 98),
      x = "Locked anchor dataset",
      y = "Genes with differential-expression support",
      caption = str_wrap("Final frozen weights are retained after reproducibility/QC steps and are not proportional to one dataset's support count.", width = 120)
    ) +
    theme_imrs_publication(base_size = 10.5) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.caption = element_text(color = "#374151", hjust = 0, margin = margin(t = 8))
    )
  save_imrs_plot(p, folder_path("Figure2_anchor_construction_weights"),
                 "Figure2A_anchor_de_summary", 10.8, 4.8, dpi = 400,
                 source_tables = required_paths$anchor_support,
                 source_code_section_or_function = "make_Figure2A",
                 notes = paste0("Restricted to the five locked anchors: ",
                                paste(LOCKED_DATASETS_MOUSE, collapse = ", "), "."))
}

prepare_weights <- function() {
  gene_col <- first_present(names(weights_tbl), c("gene", "gene_id", "ensembl_gene_id"))
  beta_col <- first_present(names(weights_tbl), c("beta_meta", "weight", "beta"))
  if (is.na(gene_col) || is.na(beta_col)) {
    stop("gene_weights.tsv must contain a gene column and a weight/beta column.", call. = FALSE)
  }
  weights <- weights_tbl %>%
    mutate(gene_id_clean = strip_ens(.data[[gene_col]]),
           applied_beta = safe_num(.data[[beta_col]])) %>%
    filter(is.finite(applied_beta))
  symbols <- symbols_tbl %>%
    transmute(gene_id_clean = strip_ens(ensembl_gene_id),
              gene_symbol = as.character(mgi_symbol)) %>%
    filter(!is.na(gene_id_clean), nzchar(gene_id_clean)) %>%
    distinct(gene_id_clean, .keep_all = TRUE)
  weights %>%
    left_join(symbols, by = "gene_id_clean") %>%
    mutate(gene_symbol = ifelse(is.na(gene_symbol) | !nzchar(gene_symbol), gene_id_clean, gene_symbol))
}

make_Figure2B <- function() {
  weights <- prepare_weights()
  retained_genes <- tibble(gene_id_clean = unique(weights$gene_id_clean))
  support_wide <- support_tbl %>%
    mutate(dataset_support_flag = logic_col(dataset_support_flag),
           gene_id_clean = strip_ens(gene_id),
           dataset_id = as.character(dataset_id)) %>%
    filter(dataset_id %in% LOCKED_DATASETS_MOUSE, gene_id_clean %in% retained_genes$gene_id_clean) %>%
    group_by(gene_id_clean, dataset_id) %>%
    summarise(supported = any(dataset_support_flag, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = dataset_id, values_from = supported, values_fill = FALSE) %>%
    right_join(retained_genes, by = "gene_id_clean")
  for (dataset_id in LOCKED_DATASETS_MOUSE) {
    if (!dataset_id %in% names(support_wide)) support_wide[[dataset_id]] <- FALSE
  }
  support_wide <- support_wide %>%
    mutate(across(all_of(LOCKED_DATASETS_MOUSE), ~ replace_na(as.logical(.x), FALSE)))
  support_matrix <- as.matrix(support_wide[, LOCKED_DATASETS_MOUSE, drop = FALSE])
  support_wide$support_n <- rowSums(support_matrix)
  support_wide$pattern_label <- apply(support_matrix, 1, function(row) {
    supported <- LOCKED_DATASETS_MOUSE[as.logical(row)]
    missing <- setdiff(LOCKED_DATASETS_MOUSE, supported)
    if (length(supported) == length(LOCKED_DATASETS_MOUSE)) {
      paste0("All ", length(LOCKED_DATASETS_MOUSE), " locked anchors")
    } else if (length(missing) == 1) {
      paste0("All except ", missing)
    } else if (length(supported) > 0) {
      paste(supported, collapse = " + ")
    } else {
      "No discovery support detected"
    }
  })
  below_threshold_n <- sum(support_wide$support_n < discovery_support_threshold, na.rm = TRUE)
  if (below_threshold_n > 0) {
    add_warning("Figure2B found ", below_threshold_n,
                " retained gene(s) below the detected discovery support threshold.")
  }
  k_locked <- length(LOCKED_DATASETS_MOUSE)
  support_label <- function(n) {
    case_when(
      n == k_locked ~ "All locked anchors",
      n == k_locked - 1L ~ "All but one locked anchor",
      n == 1L ~ "1 locked anchor",
      TRUE ~ paste0(n, " locked anchors")
    )
  }
  plot_tbl <- support_wide %>%
    mutate(
      support_n = as.integer(support_n),
      support_category = support_label(support_n)
    ) %>%
    count(support_n, support_category, name = "n_retained_genes") %>%
    filter(n_retained_genes > 0) %>%
    arrange(desc(support_n)) %>%
    mutate(support_category = factor(support_category, levels = unique(support_category)))
  pattern_notes <- support_wide %>%
    group_by(support_n, pattern_label) %>%
    summarise(n_retained_genes = n(), .groups = "drop") %>%
    mutate(
      support_category = support_label(support_n),
      supporting_datasets = purrr::map_chr(pattern_label, function(label) {
        if (startsWith(label, "All except ")) return(paste(setdiff(LOCKED_DATASETS_MOUSE, sub("^All except ", "", label)), collapse = ";"))
        if (startsWith(label, "All ")) return(paste(LOCKED_DATASETS_MOUSE, collapse = ";"))
        if (label == "No discovery support detected") return("")
        str_replace_all(label, " \\+ ", ";")
      }),
      absent_datasets = purrr::map_chr(supporting_datasets, function(supported) {
        supported_vec <- if (nzchar(supported)) unlist(strsplit(supported, ";", fixed = TRUE)) else character()
        paste(setdiff(LOCKED_DATASETS_MOUSE, supported_vec), collapse = ";")
      })
    ) %>%
    arrange(desc(support_n), pattern_label)
  readr::write_tsv(pattern_notes, figure2b_support_pattern_path, na = "NA")

  if (sum(plot_tbl$n_retained_genes, na.rm = TRUE) != nrow(retained_genes)) {
    stop("Figure2B support-count categories do not sum to retained genes.", call. = FALSE)
  }

  all_but_one_total <- plot_tbl %>%
    filter(support_n == k_locked - 1L) %>%
    summarise(total = sum(n_retained_genes), .groups = "drop") %>%
    pull(total)
  if (length(all_but_one_total) == 0L) all_but_one_total <- 0L

  missing_anchor_tbl <- pattern_notes %>%
    filter(support_n == k_locked - 1L, str_detect(pattern_label, "^All except ")) %>%
    transmute(
      missing_anchor = str_remove(pattern_label, "^All except "),
      n_retained_genes
    ) %>%
    group_by(missing_anchor) %>%
    summarise(n_retained_genes = sum(n_retained_genes), .groups = "drop") %>%
    mutate(missing_anchor = factor(missing_anchor, levels = LOCKED_DATASETS_MOUSE)) %>%
    arrange(missing_anchor)
  if (sum(missing_anchor_tbl$n_retained_genes, na.rm = TRUE) != all_but_one_total) {
    stop("Figure2B all-but-one missing-anchor breakdown does not sum to the all-but-one total.", call. = FALSE)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Figure2B combined layout requires the patchwork package.", call. = FALSE)
  }

  main_y_max <- max(plot_tbl$n_retained_genes, na.rm = TRUE)
  main_plot <- ggplot(plot_tbl, aes(x = support_category, y = n_retained_genes)) +
    geom_col(fill = "#E5E7EB", color = "#374151", linewidth = 0.45, width = 0.68) +
    geom_text(aes(label = n_retained_genes), vjust = -0.35, size = 4.2, color = "#111827") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16)), breaks = pretty_breaks(n = 5)) +
    coord_cartesian(ylim = c(0, main_y_max * 1.18), clip = "off") +
    labs(
      x = "Locked-anchor support",
      y = "Number of retained genes"
    ) +
    theme_imrs_publication(base_size = 12, legend_position = "none") +
    theme(
      axis.text.x = element_text(lineheight = 0.95, margin = margin(t = 6)),
      plot.margin = margin(10, 14, 8, 10)
    )

  mini_y_max <- if (nrow(missing_anchor_tbl) > 0L) max(missing_anchor_tbl$n_retained_genes, na.rm = TRUE) else 1
  mini_plot <- ggplot(missing_anchor_tbl, aes(x = missing_anchor, y = n_retained_genes, fill = missing_anchor)) +
    geom_col(color = "#374151", linewidth = 0.35, width = 0.7) +
    geom_text(aes(label = ifelse(n_retained_genes > 0, n_retained_genes, "")),
              vjust = -0.25, size = 3.1, color = "#111827") +
    scale_fill_manual(values = anchor_palette[LOCKED_DATASETS_MOUSE], drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)), breaks = pretty_breaks(n = 4)) +
    coord_cartesian(ylim = c(0, max(1, mini_y_max) * 1.22), clip = "off") +
    labs(
      title = "Breakdown of\nall-but-one support",
      x = "Missing anchor",
      y = "Retained genes"
    ) +
    theme_imrs_publication(base_size = 9.2, legend_position = "none") +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.8),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.title = element_text(size = 8.8, lineheight = 0.95),
      plot.margin = margin(20, 8, 8, 8)
    )

  p <- patchwork::wrap_plots(main_plot, mini_plot, nrow = 1, widths = c(2.35, 1)) +
    patchwork::plot_annotation(
      title = "Core-gene reproducibility across locked acute anchors",
      subtitle = "Retained genes show high reproducibility across the locked anchor set; all-but-one patterns are summarized separately.",
      caption = paste0(
        "Exact support-pattern counts are written to the companion support-pattern table. Locked anchors: ",
        paste(LOCKED_DATASETS_MOUSE, collapse = ", "), "."
      )
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 15, color = "#111827", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10.8, color = "#374151", lineheight = 1.08, margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, color = "#374151", hjust = 0, margin = margin(t = 8)),
      plot.margin = margin(10, 16, 10, 12)
    )

  save_imrs_plot(p, folder_path("Figure2_anchor_construction_weights"),
                 "Figure2B_core_gene_reproducibility", 9.2, 5.2, dpi = 400,
                 source_tables = c(required_paths$anchor_support, required_paths$gene_weights),
                 source_code_section_or_function = "make_Figure2B",
                 notes = paste0("Support counts use the five locked anchors; exact support-pattern counts are written to ",
                                figure2b_support_pattern_path, "."))

  detail_tbl <- pattern_notes %>%
    mutate(
      pattern_display_raw = case_when(
        support_n == k_locked ~ "All locked anchors",
        str_detect(pattern_label, "^All except ") ~ paste0("All except ", str_remove(pattern_label, "^All except ")),
        TRUE ~ pattern_label
      ),
      pattern_display = str_wrap(pattern_display_raw, width = 18)
    ) %>%
    arrange(desc(support_n), pattern_label) %>%
    mutate(pattern_display = factor(pattern_display, levels = unique(pattern_display)))
  detail_y_max <- max(detail_tbl$n_retained_genes, na.rm = TRUE)

  detail_bars <- ggplot(detail_tbl, aes(x = pattern_display, y = n_retained_genes)) +
    geom_col(fill = "#E5E7EB", color = "#374151", linewidth = 0.4, width = 0.68) +
    geom_text(aes(label = n_retained_genes), vjust = -0.3, size = 3.5, color = "#111827") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)), breaks = pretty_breaks(n = 5)) +
    coord_cartesian(ylim = c(0, detail_y_max * 1.2), clip = "off") +
    labs(x = NULL, y = "Number of retained genes") +
    theme_imrs_publication(base_size = 10.5, legend_position = "none") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(8, 12, 2, 10)
    )

  anchor_marker_tbl <- tidyr::expand_grid(
    detail_tbl %>% select(pattern_display, supporting_datasets),
    anchor_id = LOCKED_DATASETS_MOUSE
  ) %>%
    mutate(
      anchor_id = factor(anchor_id, levels = rev(LOCKED_DATASETS_MOUSE)),
      supported = map2_lgl(supporting_datasets, as.character(anchor_id), function(supported, anchor) {
        supported_vec <- if (nzchar(supported)) unlist(strsplit(supported, ";", fixed = TRUE)) else character()
        anchor %in% supported_vec
      })
    )

  detail_anchor_map <- ggplot(anchor_marker_tbl, aes(x = pattern_display, y = anchor_id)) +
    geom_tile(fill = "#F3F4F6", color = "white", linewidth = 0.45) +
    geom_tile(data = filter(anchor_marker_tbl, supported),
              aes(fill = anchor_id), color = "white", linewidth = 0.45) +
    scale_fill_manual(values = anchor_palette[LOCKED_DATASETS_MOUSE], name = "Locked anchor", drop = FALSE) +
    labs(x = "Support pattern", y = NULL, fill = "Locked anchor") +
    theme_imrs_publication(base_size = 9.5) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 8),
      legend.position = "right",
      plot.margin = margin(2, 12, 8, 10)
    )

  detail_plot <- patchwork::wrap_plots(detail_bars, detail_anchor_map,
                                       ncol = 1, heights = c(2.15, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = "Detailed locked-anchor support patterns among retained genes"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = "#111827", margin = margin(b = 6)),
      plot.margin = margin(10, 14, 10, 12)
    )

  save_imrs_plot(detail_plot, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS_anchor_support_patterns_detailed", 8.8, 5.4, dpi = 400,
                 source_tables = c(required_paths$anchor_support, required_paths$gene_weights),
                 source_code_section_or_function = "make_Figure2B",
                 notes = paste0("Detailed locked-anchor support-pattern counts are written to ",
                                figure2b_support_pattern_path, "."))
}

make_Figure2C <- function() {
  weights <- prepare_weights() %>%
    arrange(desc(abs(applied_beta))) %>%
    slice_head(n = 25) %>%
    mutate(
      gene_symbol = factor(gene_symbol, levels = rev(gene_symbol))
    )
  p <- ggplot(weights, aes(x = abs(applied_beta), y = gene_symbol)) +
    geom_col(width = 0.72, fill = "#6B7280") +
    labs(
      title = "Largest IMRS weights highlight acute discovery response genes",
      subtitle = "Genes are ranked by magnitude of the frozen weight, not signed direction.",
      x = "Absolute frozen IMRS weight",
      y = "Gene"
    ) +
    theme_imrs_publication(legend_position = "none")
  save_imrs_plot(p, folder_path("Figure2_anchor_construction_weights"),
                 "Figure2C_top_weighted_genes", 7.2, 5.2, dpi = 400,
                 source_tables = c(required_paths$gene_weights, required_paths$gene_symbols),
                 source_code_section_or_function = "make_Figure2C",
                 notes = "Panel ranks absolute frozen IMRS weight magnitude; bar fill is intentionally neutral.")
}

make_Figure2D <- function() {
  weights <- prepare_weights()
  p <- ggplot(weights, aes(x = applied_beta)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_histogram(bins = 40, fill = "#4E79A7", color = "white", linewidth = 0.2) +
    labs(
      title = "Distribution of frozen IMRS gene weights",
      subtitle = paste0("Distribution of fixed gene weights used for sample scoring; n = ", nrow(weights), " genes."),
      x = "Frozen IMRS gene weight",
      y = "Number of genes"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure2_anchor_construction_weights"),
                 "Figure2D_weight_distribution", 7, 4.6, dpi = 400,
                 source_tables = required_paths$gene_weights,
                 source_code_section_or_function = "make_Figure2D")
}

make_Figure2E <- function() {
  weights <- prepare_weights()
  count_gene_flags <- function(tbl, flag_col) {
    tbl %>%
      mutate(gene_id_clean = strip_ens(gene_id),
             is_core_gene = if ("is_core" %in% names(tbl)) logic_col(is_core) else TRUE,
             flagged = logic_col(.data[[flag_col]])) %>%
      filter(is_core_gene, !is.na(gene_id_clean), nzchar(gene_id_clean)) %>%
      group_by(gene_id_clean) %>%
      summarise(flagged = any(flagged, na.rm = TRUE), .groups = "drop") %>%
      summarise(n = sum(flagged, na.rm = TRUE), .groups = "drop") %>%
      pull(n)
  }
  heterogeneity_n <- count_gene_flags(heterogeneity_tbl, "heterogeneity_flag")
  low_power_n <- count_gene_flags(power_tbl, "low_power_flag")
  category_levels <- c(
    "Final retained scoring genes",
    "Genes with heterogeneity annotation",
    if (low_power_n > 0) "Genes with low-power annotation"
  )
  qc_palette <- c(
    "Final retained scoring genes" = "#4E79A7",
    "Genes with heterogeneity annotation" = "#F28E2B",
    "Genes with low-power annotation" = "#9CA3AF"
  )
  flag_counts <- tibble(
    category = c("Final retained scoring genes", "Genes with heterogeneity annotation",
                 if (low_power_n > 0) "Genes with low-power annotation"),
    count = c(n_distinct(weights$gene_id_clean), heterogeneity_n,
              if (low_power_n > 0) low_power_n)
  ) %>%
    mutate(category = factor(category, levels = rev(category_levels)))
  low_power_caption <- if (low_power_n == 0) {
    "No genes were removed by the low-power filter under the current criteria."
  } else {
    NA_character_
  }
  p <- ggplot(flag_counts, aes(x = count, y = category, fill = category)) +
    geom_col(width = 0.64) +
    geom_text(aes(label = count), hjust = -0.15, size = 3.4, color = "#111827") +
    scale_fill_manual(values = qc_palette, guide = "none", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = "Core-gene QC annotation summary",
      subtitle = str_wrap("QC annotations are reported for transparency and are not a mutually exclusive filtering funnel; retained genes are the final frozen scoring set from the five locked anchors.", width = 98),
      x = "Number of genes",
      y = "Gene category",
      caption = low_power_caption
    ) +
    coord_cartesian(clip = "off") +
    theme_imrs_publication(base_size = 10.5) +
    theme(
      plot.caption = element_text(color = "#374151", hjust = 0, margin = margin(t = 8)),
      plot.subtitle = element_text(color = "#374151", lineheight = 1.05, margin = margin(b = 8)),
      plot.margin = margin(8, 28, 8, 10)
    )
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3D_anchor_core_gene_qc_supplement", 8.8, 4.9, dpi = 400,
                 source_tables = c(required_paths$gene_weights, required_paths$gene_heterogeneity, required_paths$gene_power),
                 source_code_section_or_function = "make_Figure2E",
                 notes = "Moved out of main Figure 2. Standalone QC annotation summary uses distinct gene-level flag counts; categories may overlap with retained scoring genes.")
}

make_Figure3A <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% manuscript_plot_groups)
  p <- ggplot(plot_tbl, aes(x = manuscript_group, y = delta_mean_imrs_z, fill = manuscript_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_violin(alpha = 0.55, width = 0.85, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.75) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(-2, -1, 0, 1, 2, 5, 10, 25, 50)) +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "Primary validation datasets show consistent IMRS elevation",
      subtitle = "Each dot is a split contrast, not a dataset; pseudo-log scaling reduces anchor outlier compression.",
      x = "Manuscript analysis group",
      y = "Mean ΔIMRSz (pseudo-log scale)",
      fill = "Manuscript analysis group"
    ) +
    coord_cartesian(clip = "off") +
    theme_imrs_publication(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1),
      legend.box.margin = margin(t = 4),
      plot.margin = margin(8, 18, 8, 10)
    )
  save_imrs_plot(p, folder_path("Figure3_primary_validation_overview"),
                 "Figure3A_delivery_response_by_analysis_group", 8.8, 5.2, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3A",
                 notes = "Kept as Figure3A overview with pseudo-log y-axis; each dot is a split contrast.")
}

validation_tbl <- function() {
  role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% c("Primary acute validation", "Extended validation")) %>%
    mutate(manuscript_group = factor(as.character(manuscript_group),
                                     levels = c("Primary acute validation", "Extended validation")))
}

primary_validation_plot <- function(slide = FALSE) {
  tbl <- validation_tbl()
  y_range <- range(tbl$delta_mean_imrs_z, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  ann_y <- max(tbl$delta_mean_imrs_z, na.rm = TRUE) + max(0.8, 0.12 * y_span)
  ann <- tbl %>%
    group_by(manuscript_group) %>%
    summarise(n = n(),
              mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              prop_pos = mean(delta_mean_imrs_z > 0, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(label = paste0("n=", n, "\nmean=", sprintf("%.2f", mean_delta),
                          "\npositive=", scales::percent(prop_pos, accuracy = 1)),
           y = ann_y)
  ggplot(tbl, aes(x = manuscript_group, y = delta_mean_imrs_z, fill = manuscript_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.70) +
    geom_jitter(width = 0.13, height = 0, color = "#111827", size = 2, alpha = 0.75, show.legend = FALSE) +
    geom_text(data = ann, aes(x = manuscript_group, y = y, label = label),
              inherit.aes = FALSE, size = ifelse(slide, 4.1, 3.8),
              fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
    labs(
      title = ifelse(slide, "Primary and extended datasets show positive IMRS shifts", "Primary acute validation supports delivery-associated IMRS elevation"),
      subtitle = ifelse(slide, "Each point is a split contrast grouped for slide-ready comparison.", "Each point is a split contrast; the dashed line marks no delivery-associated score shift."),
      x = "Manuscript analysis group",
      y = "Mean delivery-minus-control IMRS z-score",
      fill = "Manuscript analysis group"
    ) +
    coord_cartesian(clip = "off") +
    theme_imrs_publication(base_size = ifelse(slide, 12.5, 11.5), legend_position = "none") +
    theme(plot.margin = margin(8, 24, 8, 10))
}

make_Figure3B <- function() {
  save_imrs_plot(primary_validation_plot(FALSE), folder_path("Figure3_primary_validation_overview"),
                 "Figure3B_primary_validation_summary", 7.2, 4.8, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3B")
}

make_Figure3C <- function() {
  save_imrs_plot(primary_validation_plot(TRUE), folder_path("Figure3_primary_validation_overview"),
                 "Figure3C_primary_extended_slide_summary", 8.5, 5.5, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3C",
                 notes = "Canonical larger slide version of Figure3B.")
}

make_Figure3D <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% manuscript_plot_groups) %>%
    mutate(
      contrast_short = short_text(display_text(split_id), 26),
      tissue_short = case_when(
        str_to_lower(tissue) == "draining lymph nodes" ~ "dLN (draining lymph node)",
        str_to_lower(tissue) == "dln" ~ "dLN (draining lymph node)",
        TRUE ~ display_text(tissue)
      ),
      label = short_text(paste(dataset_id, tissue_short, format_time_label(time_h), contrast_short, sep = " | "), 72)
    ) %>%
    arrange(manuscript_group, delta_mean_imrs_z) %>%
    mutate(label = factor(label, levels = unique(as.character(label))))
  p <- ggplot(plot_tbl, aes(x = delta_mean_imrs_z, y = label, color = manuscript_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_point(size = 1.35, alpha = 0.9) +
    facet_grid(manuscript_group ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Contrast-level validation responses are directionally positive",
      subtitle = "Full contrast forest; dLN means draining lymph node.",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "Contrast label",
      color = "Manuscript analysis group"
    ) +
    theme_imrs_publication(base_size = 6.8, legend_position = "none") +
    theme(strip.text.y = element_text(angle = 0))
  save_imrs_plot(p, folder_path("FigureS2_detailed_validation_forests"),
                 "FigureS2A_full_validation_contrast_forest", 7.8, max(8.8, 0.16 * nrow(plot_tbl) + 2),
                 dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3D",
                 notes = "Dense appendix version retained; see Figure3D_simplified_top_contrasts for slide-friendly view.")
}

make_Figure3D_simplified <- function() {
  base_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% manuscript_plot_groups)
  locked_anchor_tbl <- base_tbl %>%
    filter(
      as.character(manuscript_group) == "Locked anchor",
      dataset_id %in% LOCKED_DATASETS_MOUSE,
      is.na(time_h) | time_h <= 24
    ) %>%
    group_by(dataset_id) %>%
    slice_max(order_by = abs(delta_mean_imrs_z), n = 1, with_ties = FALSE) %>%
    ungroup()
  other_tbl <- base_tbl %>%
    filter(as.character(manuscript_group) != "Locked anchor") %>%
    group_by(manuscript_group) %>%
    slice_max(order_by = abs(delta_mean_imrs_z), n = 4, with_ties = FALSE) %>%
    ungroup()
  plot_tbl <- bind_rows(locked_anchor_tbl, other_tbl) %>%
    mutate(
      label = short_text(dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE), 48),
      label = ordered_factor(label, delta_mean_imrs_z)
    )
  locked_check <- tibble(dataset_id = LOCKED_DATASETS_MOUSE, expected_locked_anchor = TRUE) %>%
    left_join(
      plot_tbl %>%
        filter(as.character(manuscript_group) == "Locked anchor",
               dataset_id %in% LOCKED_DATASETS_MOUSE) %>%
        transmute(
          dataset_id,
          appears_in_Figure4A = TRUE,
          plotted_group = as.character(manuscript_group),
          plotted_label = as.character(label),
          mean_delta = delta_mean_imrs_z
        ),
      by = "dataset_id"
    ) %>%
    mutate(
      appears_in_Figure4A = replace_na(appears_in_Figure4A, FALSE),
      plotted_group = ifelse(appears_in_Figure4A, plotted_group, NA_character_),
      plotted_label = ifelse(appears_in_Figure4A, plotted_label, NA_character_),
      reason_if_missing = ifelse(appears_in_Figure4A, "included",
                                 "missing from locked-anchor Figure4A panel")
    ) %>%
    arrange(match(dataset_id, LOCKED_DATASETS_MOUSE))
  readr::write_tsv(locked_check, figure4a_anchor_inclusion_check_path, na = "NA")
  if (any(!locked_check$appears_in_Figure4A)) {
    stop("Figure4A locked-anchor inclusion check failed; missing: ",
         paste(locked_check$dataset_id[!locked_check$appears_in_Figure4A], collapse = ", "),
         call. = FALSE)
  }
  p <- ggplot(plot_tbl, aes(x = delta_mean_imrs_z, y = label, color = manuscript_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_point(size = 2.6, alpha = 0.95) +
    facet_wrap(~ manuscript_group, scales = "free_y", ncol = 2) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(
      title = "Validation and anchor responses remain positive across groups",
      subtitle = "Five locked anchors appear once; validation panels show the largest score shifts per group.",
      x = "Mean delivery-minus-control ΔIMRSz",
      y = "Dataset context",
      color = "Manuscript analysis group"
    ) +
    coord_cartesian(clip = "off") +
    theme_imrs_publication(base_size = 10.5, legend_position = "none") +
    theme(
      axis.text.y = element_text(size = 9.5),
      strip.text = element_text(size = 11),
      plot.margin = margin(8, 22, 8, 10)
    )
  save_imrs_plot(p, folder_path("Figure4_validation_detail_and_discrimination"),
                 "Figure4A_top_contrast_responses", 10, 6.8, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3D_simplified",
                 notes = paste0("Slide-friendly summary with all five locked anchors represented once; inclusion check: ",
                                figure4a_anchor_inclusion_check_path))
}

make_Figure3E <- function() {
  plot_tbl <- role_pass_for_plot %>% filter(is.finite(auc_imrs_z))
  p <- ggplot(plot_tbl, aes(x = auc_imrs_z, fill = manuscript_group)) +
    geom_histogram(bins = 20, alpha = 0.78, color = "white", position = "identity") +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "Sample-level AUC provides secondary validation evidence",
      subtitle = "Area under the ROC curve (AUC) is secondary evidence; delivery-minus-control IMRS shift remains the primary metric.",
      x = "AUC for IMRS z-score separation",
      y = "Number of split contrasts",
      fill = "Manuscript analysis group"
    ) +
    theme_imrs_publication(base_size = 10.5) +
    theme(plot.margin = margin(8, 18, 8, 10))
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3E_auc_secondary_validation", 8.3, 4.8, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3E",
                 notes = "Moved from main validation-detail figure to supplement because AUC is secondary evidence.")
}

make_Figure3F <- function() {
  plot_tbl <- role_pass_for_plot %>% filter(is.finite(delta_mean_imrs_z), is.finite(auc_imrs_z))
  p <- ggplot(plot_tbl, aes(x = delta_mean_imrs_z, y = auc_imrs_z,
                            color = manuscript_group, size = n_controls + n_delivery)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_hline(yintercept = 0.5, linetype = "dotted", linewidth = 0.35, color = "#4B5563") +
    geom_jitter(width = 0.05, height = 0.008, alpha = 0.58) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    scale_size_continuous(range = c(1.5, 5)) +
    scale_y_continuous(breaks = c(0.5, 0.75, 0.9, 1.0)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "Contrast-level IMRS shifts align with sample-level separation",
      subtitle = "AUC is secondary evidence, not the primary IMRS endpoint; small jitter reveals overlapping contrasts near AUC = 1.",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "AUC for control–delivery separation",
      color = "Manuscript analysis group",
      size = "Sample size"
    ) +
    coord_cartesian(ylim = c(0.35, 1.03)) +
    theme_imrs_publication(base_size = 10.5) +
    theme(plot.margin = margin(8, 24, 8, 10))
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3F_delta_auc_secondary_relationship", 8.4, 5.1, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_Figure3F",
                 notes = "Moved from main validation-detail figure to supplement; jitter added to show overlap near AUC = 1.")
}

make_Figure4A <- function() {
  plot_tbl <- perm_summary_tbl %>%
    filter(permutation_status == "ok") %>%
    attach_role_group() %>%
    mutate(order_index = row_number()) %>%
    arrange(observed_delta_mean_imrs_z) %>%
    mutate(order_index = row_number())
  p <- ggplot(plot_tbl, aes(x = order_index)) +
    geom_linerange(aes(ymin = null_q025_delta, ymax = null_q975_delta),
                   color = "#BDBDBD", linewidth = 0.7) +
    geom_point(aes(y = observed_delta_mean_imrs_z, color = manuscript_group),
               size = 1.7, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Observed IMRS shifts exceed label-permutation expectations",
      subtitle = "Observed delivery-minus-control shifts are compared with 95% within-contrast label-permutation intervals.",
      x = "Split contrasts ordered by observed delivery-minus-control IMRS z-score",
      y = "Mean delivery-minus-control IMRS z-score",
      color = "Manuscript analysis group"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure5_permutation_null_analysis"),
                 "Figure5A_label_permutation_observed_vs_null", 8.4, 5.1, dpi = 400,
                 source_tables = required_paths$label_permutation_summary,
                 source_code_section_or_function = "make_Figure4A")
}

make_Figure4B <- function() {
  plot_tbl <- perm_summary_tbl %>%
    filter(permutation_status == "ok") %>%
    attach_role_group()
  p <- ggplot(plot_tbl, aes(x = empirical_p_greater)) +
    geom_histogram(binwidth = 0.05, boundary = 0, color = "white", fill = "#4E79A7") +
    facet_wrap(~ manuscript_group, ncol = 2) +
    labs(
      title = "Permutation p-values support delivery-associated IMRS elevation",
      subtitle = "Empirical greater-than p-values test whether delivery samples have higher control-standardized IMRS z-scores than controls.",
      x = "Empirical p-value for greater-than null",
      y = "Number of split contrasts"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure5_permutation_null_analysis"),
                 "Figure5B_empirical_pvalue_distribution", 8, 5.4, dpi = 400,
                 source_tables = required_paths$label_permutation_summary,
                 source_code_section_or_function = "make_Figure4B")
}

make_Figure4C <- function() {
  plot_tbl <- perm_summary_tbl %>%
    filter(permutation_status == "ok") %>%
    attach_role_group()
  p <- ggplot(plot_tbl, aes(x = manuscript_group, y = observed_delta_mean_imrs_z, fill = manuscript_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72) +
    geom_jitter(width = 0.13, height = 0, size = 1.7, alpha = 0.75) +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Permutation-tested IMRS shifts differ by analysis group",
      subtitle = "Primary acute validation is summarized separately from extended validation and secondary support.",
      x = "Manuscript analysis group",
      y = "Observed mean delivery-minus-control IMRS z-score",
      fill = "Manuscript analysis group"
    ) +
    theme_imrs_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_imrs_plot(p, folder_path("Figure5_permutation_null_analysis"),
                 "Figure5C_permutation_response_by_analysis_group", 8.2, 5.2, dpi = 400,
                 source_tables = required_paths$label_permutation_summary,
                 source_code_section_or_function = "make_Figure4C")
}

make_Figure4D <- function() {
  plot_tbl <- baseline_long_tbl %>%
    filter(logic_col(pass)) %>%
    attach_role_group() %>%
    mutate(score_display = map_score_label(score_label, score_id))
  p <- ggplot(plot_tbl, aes(x = manuscript_group, y = safe_num(delta_score), fill = manuscript_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72) +
    geom_jitter(width = 0.12, height = 0, size = 1.35, alpha = 0.62) +
    facet_wrap(~ score_display, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Baseline signatures provide comparator response profiles",
      subtitle = "Benchmark signatures are positive-control comparators, not replacements for IMRS.",
      x = "Manuscript analysis group",
      y = "Mean delivery-minus-control signature z-score",
      fill = "Manuscript analysis group"
    ) +
    theme_imrs_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_imrs_plot(p, folder_path("Figure6_baseline_signature_benchmarking"),
                 "Figure6A_baseline_delta_by_analysis_group", 9, 6.4, dpi = 400,
                 source_tables = required_paths$baseline_contrast_long,
                 source_code_section_or_function = "make_Figure4D")
}

make_Figure4E <- function() {
  plot_tbl <- baseline_paired_tbl %>%
    mutate(score_display = map_score_label(signature_label, signature_id),
           baseline_delta_score = safe_num(baseline_delta_score),
           imrs_delta_score = safe_num(imrs_delta_score))
  p <- ggplot(plot_tbl, aes(x = imrs_delta_score, y = baseline_delta_score, color = score_display)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_point(alpha = 0.72, size = 1.8) +
    facet_wrap(~ score_display, scales = "free_y") +
    scale_color_manual(values = score_palette, drop = FALSE) +
    labs(
      title = "IMRS and baseline signatures are compared within matched contrasts",
      subtitle = "Each point is the same split contrast, supporting benchmark interpretation without implying categorical superiority.",
      x = "IMRS delivery-minus-control z-score",
      y = "Baseline signature delivery-minus-control z-score",
      color = "Benchmark score"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure6_baseline_signature_benchmarking"),
                 "Figure6B_paired_imrs_baseline_contrast_comparison", 8.5, 5.5, dpi = 400,
                 source_tables = required_paths$baseline_paired,
                 source_code_section_or_function = "make_Figure4E")
}

make_Figure4F <- function() {
  plot_tbl <- baseline_long_tbl %>%
    filter(logic_col(pass)) %>%
    attach_role_group() %>%
    mutate(score_display = map_score_label(score_label, score_id),
           delta_score = safe_num(delta_score)) %>%
    filter(is.finite(delta_score)) %>%
    group_by(score_id, score_label, score_display, manuscript_group) %>%
    summarise(
      n_contrasts = n(),
      proportion_positive_delta = mean(delta_score > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(proportion_positive_delta))
  p <- ggplot(plot_tbl, aes(x = manuscript_group, y = proportion_positive_delta, fill = score_display)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_fill_manual(values = score_palette, drop = FALSE) +
    labs(
      title = "IMRS and benchmark signatures are directionally compared",
      subtitle = "Bars show the fraction of split contrasts with higher scores in delivery samples than controls.",
      x = "Manuscript analysis group",
      y = "Proportion of positive delivery-associated contrasts",
      fill = "Benchmark score"
    ) +
    theme_imrs_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_imrs_plot(p, folder_path("Figure6_baseline_signature_benchmarking"),
                 "Figure6C_benchmark_directionality_summary", 8.8, 5.2, dpi = 400,
                 source_tables = required_paths$baseline_contrast_long,
                 source_code_section_or_function = "make_Figure4F")
}

make_Figure4G <- function() {
  sample_role <- step09_sample_tbl %>%
    left_join(role_lookup %>% select(gse_id, split_id, manuscript_group), by = c("gse_id", "split_id")) %>%
    mutate(manuscript_group = factor(as.character(manuscript_group), levels = manuscript_group_order)) %>%
    clean_dataset_role_for_plot("manuscript_group")
  plot_tbl <- sample_role %>%
    inner_join(baseline_scores_tbl, by = c("gse_id" = "dataset_id", "sample_id"),
               relationship = "many-to-many") %>%
    mutate(score_display = map_score_label(signature_label, signature_id),
           signature_z = safe_num(signature_z),
           imrs_z = safe_num(imrs_z)) %>%
    filter(is.finite(signature_z), is.finite(imrs_z))
  p <- ggplot(plot_tbl, aes(x = signature_z, y = imrs_z, color = manuscript_group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#111827", linewidth = 0.45) +
    facet_wrap(~ score_display, scales = "free") +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Sample-level IMRS variation is benchmarked against baseline signatures",
      subtitle = "Each point is a scored sample within a split-level contrast.",
      x = "Control-standardized baseline signature z-score",
      y = "Control-standardized IMRS z-score",
      color = "Manuscript analysis group"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure6_baseline_signature_benchmarking"),
                 "Figure6D_sample_level_baseline_correlation", 8.6, 5.6, dpi = 400,
                 source_tables = c(required_paths$baseline_scores_sample, required_paths$step09_sample),
                 source_code_section_or_function = "make_Figure4G")
}

join_role_group <- function(df) {
  attach_role_group(df)
}

make_Figure5A <- function() {
  plot_tbl <- loo_tbl %>%
    join_role_group() %>%
    mutate(original_delta_mean_imrs_z = safe_num(original_delta_mean_imrs_z),
           leave_one_gene_out_delta_mean_imrs_z = safe_num(leave_one_gene_out_delta_mean_imrs_z))
  p <- ggplot(plot_tbl, aes(x = original_delta_mean_imrs_z, y = leave_one_gene_out_delta_mean_imrs_z)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "#4B5563") +
    geom_point(aes(color = manuscript_group), alpha = 0.45, size = 1.4) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "IMRS scores remain stable after single-gene removal",
      subtitle = "Each point compares the original response with the response after removing one weighted gene.",
      x = "Original mean delivery-minus-control IMRS z-score",
      y = "After single-gene removal IMRS z-score",
      color = "Analysis group"
    ) +
    theme_imrs_publication(base_size = 10.5) +
    theme(plot.margin = margin(8, 18, 8, 10))
  save_imrs_plot(p, folder_path("Figure7_gene_threshold_robustness"),
                 "Figure7A_leave_one_gene_out_delta_correlation", 8.8, 5.2, dpi = 400,
                 source_tables = required_paths$leave_one_gene_out,
                 source_code_section_or_function = "make_Figure5A")
}

make_Figure5B <- function() {
  plot_tbl <- loo_tbl %>%
    mutate(absolute_change_delta = safe_num(absolute_change_delta),
           direction_flag = factor(ifelse(logic_col(direction_preserved), "Direction preserved", "Any direction change"),
                                   levels = names(preservation_palette))) %>%
    group_by(removed_gene_id, removed_gene_symbol, removal_rank) %>%
    summarise(max_absolute_change_delta = max(absolute_change_delta, na.rm = TRUE),
              direction_flag = ifelse(any(direction_flag == "Any direction change", na.rm = TRUE),
                                      "Any direction change", "Direction preserved"),
              .groups = "drop") %>%
    arrange(desc(max_absolute_change_delta)) %>%
    slice_head(n = 20) %>%
    mutate(label = ifelse(!is.na(removed_gene_symbol) & nzchar(removed_gene_symbol),
                          removed_gene_symbol, removed_gene_id),
           label = factor(label, levels = rev(unique(label))),
           direction_flag = factor(direction_flag, levels = names(preservation_palette)))
  p <- ggplot(plot_tbl, aes(x = max_absolute_change_delta, y = label, fill = direction_flag)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = preservation_palette, drop = FALSE) +
    labs(
      title = "No single gene dominates IMRS contrast-level responses",
      subtitle = "Bars show the largest absolute response change observed after removing each weighted gene.",
      x = "Maximum absolute change in delivery-minus-control IMRS z-score",
      y = "Removed gene",
      fill = "Direction status"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure7_gene_threshold_robustness"),
                 "Figure7B_top_gene_removal_effect", 7.5, 5.5, dpi = 400,
                 source_tables = required_paths$leave_one_gene_out,
                 source_code_section_or_function = "make_Figure5B")
}

make_Figure5C <- function() {
  plot_tbl <- dominance_tbl %>%
    join_role_group() %>%
    mutate(median_max_contribution_fraction = safe_num(median_max_contribution_fraction))
  p <- ggplot(plot_tbl, aes(x = median_max_contribution_fraction, fill = manuscript_group)) +
    geom_histogram(binwidth = 0.025, alpha = 0.78, position = "identity") +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "IMRS response is not driven by a single dominant gene",
      subtitle = "Top-gene contribution fractions summarize score concentration within each contrast.",
      x = "Median top-gene contribution fraction",
      y = "Number of contrasts",
      fill = "Analysis group"
    ) +
    theme_imrs_publication(base_size = 10.5) +
    theme(plot.margin = margin(8, 18, 8, 10))
  save_imrs_plot(p, folder_path("Figure7_gene_threshold_robustness"),
                 "Figure7C_gene_dominance_distribution", 8.6, 5.1, dpi = 400,
                 source_tables = required_paths$gene_dominance,
                 source_code_section_or_function = "make_Figure5C")
}

make_Figure5D <- function() {
  plot_tbl <- threshold_summary_tbl %>%
    mutate(across(c(min_abs_log2FC, min_up_anchor_count, fdr_support,
                    external_mean_delta_imrs_z), safe_num))
  p <- ggplot(plot_tbl, aes(x = factor(min_abs_log2FC), y = factor(fdr_support),
                            fill = external_mean_delta_imrs_z)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = paste0("n=", n_genes_selected, "\n",
                                 sprintf("%.2f", external_mean_delta_imrs_z))), size = 3) +
    facet_wrap(~ min_up_anchor_count, labeller = labeller(min_up_anchor_count = label_both)) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, na.value = "grey85") +
    labs(
      title = "IMRS response is robust to gene-selection threshold variation",
      subtitle = "Mean delivery-associated score shifts remain directionally consistent across core-gene selection settings.",
      x = "Minimum absolute log2 fold change",
      y = "False discovery rate cutoff",
      fill = "Mean delivery-minus-control\nIMRS z-score"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure7_gene_threshold_robustness"),
                 "Figure7D_threshold_sensitivity_heatmap", 8.5, 5.5, dpi = 400,
                 source_tables = required_paths$threshold_summary,
                 source_code_section_or_function = "make_Figure5D")
}

make_Figure5E <- function() {
  plot_tbl <- threshold_detail_tbl %>%
    mutate(across(c(min_abs_log2FC, min_up_anchor_count, fdr_support, delta_mean_imrs_z), safe_num)) %>%
    filter(sensitivity_scope == "external_full_3_anchor_weights") %>%
    mutate(fdr_label = factor(as.character(fdr_support), levels = names(threshold_fdr_palette)))
  p <- ggplot(plot_tbl, aes(x = factor(min_abs_log2FC), y = delta_mean_imrs_z, fill = fdr_label)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
    facet_wrap(~ min_up_anchor_count, labeller = labeller(min_up_anchor_count = label_both)) +
    scale_fill_manual(values = threshold_fdr_palette, drop = FALSE) +
    labs(
      title = "Contrast-level IMRS shifts remain positive across thresholds",
      subtitle = "Distributions summarize delivery-minus-control responses across validation contrasts for each selection setting.",
      x = "Minimum absolute log2 fold change",
      y = "Mean delivery-minus-control IMRS z-score",
      fill = "FDR cutoff"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure7_gene_threshold_robustness"),
                 "Figure7E_threshold_sensitivity_delta_distribution", 8.5, 5.5, dpi = 400,
                 source_tables = required_paths$threshold_detail,
                 source_code_section_or_function = "make_Figure5E")
}

make_Figure5F <- function() {
  source_audit_one <- function(path) {
    file_exists <- file.exists(path)
    if (!file_exists) {
      return(tibble(
        source_file = path,
        file_exists = FALSE,
        n_rows = NA_integer_,
        holdout_anchor_values = NA_character_,
        n_unique_holdout_anchors = NA_integer_,
        contains_GSE39129 = FALSE,
        contains_GSE167521 = FALSE,
        contains_GSE264344 = FALSE,
        contains_GSE279372 = FALSE,
        contains_GSE279744 = FALSE,
        supports_five_anchor_LOAO = FALSE,
        notes = "file missing"
      ))
    }
    df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, guess_max = 100000)
    holdout_values <- if ("holdout_anchor" %in% names(df)) {
      sort(unique(as.character(df$holdout_anchor)))
    } else {
      character()
    }
    holdout_values <- holdout_values[!is.na(holdout_values) & nzchar(holdout_values)]
    contains <- LOCKED_DATASETS_MOUSE %in% holdout_values
    tibble(
      source_file = path,
      file_exists = TRUE,
      n_rows = nrow(df),
      holdout_anchor_values = paste(holdout_values, collapse = ";"),
      n_unique_holdout_anchors = length(holdout_values),
      contains_GSE39129 = contains[[1]],
      contains_GSE167521 = contains[[2]],
      contains_GSE264344 = contains[[3]],
      contains_GSE279372 = contains[[4]],
      contains_GSE279744 = contains[[5]],
      supports_five_anchor_LOAO = setequal(holdout_values, LOCKED_DATASETS_MOUSE),
      notes = case_when(
        !"holdout_anchor" %in% names(df) ~ "no holdout_anchor column",
        setequal(holdout_values, LOCKED_DATASETS_MOUSE) ~ "contains all five locked production anchors",
        TRUE ~ "does not contain all five locked production anchors"
      )
    )
  }
  old_loao_detail <- file.path(audit_results_dir, "leave_one_anchor_out_contrast_details.tsv")
  old_loao_summary <- file.path(audit_results_dir, "leave_one_anchor_out_summary.tsv")
  loao_source_audit <- bind_rows(
    source_audit_one(old_loao_detail),
    source_audit_one(old_loao_summary),
    source_audit_one(required_paths$leave_one_anchor_detail),
    source_audit_one(required_paths$leave_one_anchor_summary)
  )
  readr::write_tsv(loao_source_audit, figure8a_loao_source_audit_path, na = "NA")
  holdout_values <- loao_detail_tbl %>%
    pull(holdout_anchor) %>%
    as.character() %>%
    unique() %>%
    sort()
  can_support_five_anchor <- setequal(holdout_values, locked_anchor_ids)
  if (!can_support_five_anchor) {
    missing_holdouts <- setdiff(locked_anchor_ids, holdout_values)
    stop("Figure8A five-anchor LOAO source is missing holdout anchor(s): ",
         paste(missing_holdouts, collapse = ", "),
         ". Source: ", required_paths$leave_one_anchor_detail, call. = FALSE)
  }
  holdout_levels <- locked_anchor_ids
  plot_tbl <- loao_detail_tbl %>%
    mutate(delta_mean_imrs_z = safe_num(delta_mean_imrs_z),
           holdout_anchor = factor(holdout_anchor, levels = holdout_levels))
  p <- ggplot(plot_tbl, aes(x = holdout_anchor, y = delta_mean_imrs_z, color = holdout_anchor)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_jitter(width = 0.16, height = 0, size = 2, alpha = 0.75) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "#111827", linewidth = 0.35) +
    scale_color_manual(values = anchor_palette[holdout_levels], breaks = holdout_levels, drop = FALSE) +
    labs(
      title = "Five-anchor leave-one-anchor-out analysis preserves positive IMRS responses",
      subtitle = "Each point shows a held-out anchor contrast scored with weights reconstructed from the other four locked anchors.",
      x = "Held-out locked anchor dataset",
      y = "Mean delivery-minus-control IMRS z-score",
      color = "Held-out anchor"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure8_anchor_coefficient_robustness"),
                 "Figure8A_leave_one_anchor_out_delta", 7.2, 5, dpi = 400,
                 source_tables = required_paths$leave_one_anchor_detail,
                 source_code_section_or_function = "make_Figure5F",
                 notes = paste0("Source table supports true five-anchor leave-one-anchor-out plotting. Source audit: ",
                                figure8a_loao_source_audit_path))
}

make_Figure5G <- function() {
  add_warning("Figure5G sample-level leave-one-anchor-out scores were not present as a result table; using existing Step09 anchor sample-level IMRS z-scores for the sample-distribution companion panel.")
  sample_check <- step09_sample_tbl %>%
    filter(gse_id %in% locked_anchor_ids) %>%
    mutate(condition_simple = as.character(condition_simple)) %>%
    group_by(gse_id) %>%
    summarise(
      n_sample_rows = n(),
      n_control = sum(condition_simple == "CONTROL", na.rm = TRUE),
      n_delivery = sum(condition_simple == "DELIVERY", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(tibble(gse_id = locked_anchor_ids), by = "gse_id") %>%
    mutate(
      n_sample_rows = replace_na(n_sample_rows, 0L),
      n_control = replace_na(n_control, 0L),
      n_delivery = replace_na(n_delivery, 0L),
      included_in_plot = n_sample_rows > 0 & n_control > 0 & n_delivery > 0,
      reason_if_missing = case_when(
        included_in_plot ~ "included",
        n_sample_rows == 0 ~ "no sample-level rows in step09_split_sample_level.tsv",
        n_control == 0 ~ "no CONTROL sample rows",
        n_delivery == 0 ~ "no DELIVERY sample rows",
        TRUE ~ "not included"
      )
    ) %>%
    arrange(match(gse_id, locked_anchor_ids))
  readr::write_tsv(sample_check, figure8b_sample_source_check_path, na = "NA")
  plot_tbl <- step09_sample_tbl %>%
    filter(gse_id %in% locked_anchor_ids, condition_simple %in% names(condition_palette)) %>%
    left_join(step09_eval_tbl %>% select(gse_id, split_id, tissue, time_h) %>% distinct(),
              by = c("gse_id", "split_id")) %>%
    mutate(
      holdout_anchor = factor(gse_id, levels = locked_anchor_ids),
      condition_simple = factor(condition_simple, levels = names(condition_palette)),
      imrs_z = safe_num(imrs_z)
    ) %>%
    filter(is.finite(imrs_z))
  p <- ggplot(plot_tbl, aes(x = condition_simple, y = imrs_z, fill = condition_simple)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted", color = "#4B5563") +
    geom_violin(alpha = 0.55, trim = FALSE, width = 0.8) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(aes(color = holdout_anchor), width = 0.08, size = 0.8, alpha = 0.45) +
    facet_wrap(~ holdout_anchor, scales = "free_y") +
    scale_x_discrete(labels = condition_labels) +
    scale_fill_manual(values = condition_palette, labels = condition_labels, drop = FALSE) +
    scale_color_manual(values = anchor_palette[locked_anchor_ids], breaks = locked_anchor_ids,
                       drop = FALSE, guide = "none") +
    labs(
      title = "Locked anchor datasets retain delivery-associated sample-level separation",
      subtitle = "Control-standardized IMRS z-scores are shown across the five production locked anchor datasets.",
      x = "Sample condition",
      y = "control-standardized IMRS z-score",
      fill = "Sample condition"
    ) +
    theme_imrs_publication()
  save_imrs_plot(p, folder_path("Figure8_anchor_coefficient_robustness"),
                 "Figure8B_anchor_sample_distributions", 8.5, 5.5, dpi = 400,
                 source_tables = required_paths$step09_sample,
                 source_code_section_or_function = "make_Figure5G",
                 notes = paste0("No leave-one-anchor-out sample-level result table was available; this panel is regenerated from existing Step09 anchor sample-level scores and uses standardized CONTROL/DELIVERY colors. Source check: ",
                                figure8b_sample_source_check_path))
}

make_Figure5H <- function() {
  plot_tbl <- dominance_tbl %>%
    join_role_group() %>%
    mutate(median_max_contribution_fraction = safe_num(median_max_contribution_fraction))
  p <- ggplot(plot_tbl, aes(x = median_max_contribution_fraction, fill = manuscript_group)) +
    geom_histogram(binwidth = 0.025, alpha = 0.78, position = "identity") +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Top-gene contribution remains low across contrasts",
      subtitle = "Most contrasts have low median top-gene contribution, supporting a non-degenerate score.",
      x = "Median top-gene contribution fraction",
      y = "Number of split contrasts",
      fill = "Manuscript analysis group"
    ) +
    theme_imrs_publication(base_size = 10)
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3C_appendix_gene_dominance", 7, 3.8, dpi = 400,
                 source_tables = required_paths$gene_dominance,
                 source_code_section_or_function = "make_Figure5H")
}

make_Figure5I <- function() {
  plot_tbl <- coefficient_summary_tbl %>%
    filter(metric %in% c(
      "leave_one_gene_out_direction_preserved_fraction",
      "threshold_sensitivity_external_positive_delta_min",
      "gene_dominance_median_max_contribution_fraction",
      "gene_dominance_max_max_contribution_fraction"
    )) %>%
    mutate(
      value_num = safe_num(value),
      metric_label = recode(
        metric,
        leave_one_gene_out_direction_preserved_fraction =
          "Direction preserved after single-gene removal",
        threshold_sensitivity_external_positive_delta_min =
          "Minimum positive-response fraction across thresholds",
        gene_dominance_median_max_contribution_fraction =
          "Median top-gene contribution fraction",
        gene_dominance_max_max_contribution_fraction =
          "Largest top-gene contribution fraction"
      )
    ) %>%
    filter(is.finite(value_num))
  p <- ggplot(plot_tbl, aes(x = reorder(metric_label, value_num), y = value_num)) +
    geom_col(width = 0.65, fill = "#4E79A7") +
    coord_flip() +
    labs(
      title = "Frozen IMRS coefficient robustness summary",
      subtitle = "Single-gene removal, threshold variation, and gene-dominance metrics.",
      x = "Robustness and sensitivity metric",
      y = "Sensitivity metric value"
    ) +
    theme_imrs_publication(base_size = 10.5)
  save_imrs_plot(p, folder_path("Figure8_anchor_coefficient_robustness"),
                 "Figure8C_coefficient_sensitivity_summary", 8.8, 4.8, dpi = 400,
                 source_tables = required_paths$coefficient_sensitivity_summary,
                 source_code_section_or_function = "make_Figure5I")
}

make_FigureSA <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% validation_plot_groups) %>%
    mutate(validation_group = factor(as.character(manuscript_group),
                                     levels = validation_plot_groups))
  p <- ggplot(plot_tbl, aes(x = validation_group, y = delta_mean_imrs_z, fill = validation_group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.58) +
    geom_jitter(width = 0.14, height = 0, size = 2.1, alpha = 0.85) +
    scale_fill_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "Acute validation shows stronger delivery-associated IMRS elevation",
      subtitle = "Contrast-level score shifts are compared across validation timing groups.",
      x = "Validation timing group",
      y = "Mean delivery-minus-control IMRS z-score",
      fill = "Validation group"
    ) +
    theme_imrs_publication(legend_position = "none") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1A_acute_vs_late_delta_boxplot", 7.5, 5, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSA")
}

make_FigureSB <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) != "Locked anchor") %>%
    arrange(manuscript_group, dataset_id, delta_mean_imrs_z) %>%
    mutate(
      contrast_short = short_text(paste(dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE),
                                        display_text(split_id), sep = " | "), 78),
      contrast_short = factor(contrast_short, levels = unique(as.character(contrast_short)))
    )
  p <- ggplot(plot_tbl, aes(x = delta_mean_imrs_z, y = contrast_short, color = manuscript_group)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_point(size = 1.7, alpha = 0.9) +
    facet_grid(manuscript_group ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    labs(
      title = "External validation responses vary by timing and biological context",
      subtitle = "Full forest plot; dLN means draining lymph node. See simplified by-dataset summary for slides.",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "Contrast",
      color = "Validation group"
    ) +
    theme_imrs_publication(base_size = 7, legend_position = "none") +
    theme(strip.text.y = element_text(angle = 0))
  save_imrs_plot(p, folder_path("FigureS2_detailed_validation_forests"),
                 "FigureS2B_acute_vs_late_full_forest", 9, max(6, 0.18 * nrow(plot_tbl) + 2), dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSB",
                 notes = "Dense appendix version retained; see FigureSB_simplified_by_dataset.")
}

make_FigureSB_simplified <- function() {
  corrected_split_id <- paste0(
    "GSE262515_design__T=cell_line__H=16__B=NA__",
    "G=delivery_tlnp_sinc_50nm__VS=baseline_0"
  )
  corrected_row <- weak_tbl %>%
    filter(
      dataset_id == "GSE262515_cell_line",
      treatment_group == "delivery_tlnp_sinc_50nm",
      split_id == corrected_split_id
    )
  if (nrow(corrected_row) != 1L) {
    stop("Supplementary Figure S1B requires exactly one matching GSE262515 siNC cell-line row.",
         call. = FALSE)
  }

  category_levels <- c(
    "Disease-rescue model",
    "Distal/adaptive tissue context",
    "Late timepoint",
    "Low-inflammatory formulation design",
    "Therapeutic cargo/context effect",
    "Tissue/time kinetic effect"
  )
  support_levels <- c("low", "medium", "high")
  support_palette <- c(low = "#4E79A7", medium = "#F28E2B", high = "#E15759")
  corrected_flag <- weak_tbl$dataset_id == "GSE262515_cell_line" &
    weak_tbl$treatment_group == "delivery_tlnp_sinc_50nm" &
    weak_tbl$split_id == corrected_split_id

  row_tbl <- weak_tbl %>%
    mutate(
      explanation_category_plot = if_else(
        corrected_flag,
        "therapeutic_cargo_specific_effect",
        as.character(explanation_category)
      ),
      support_level = if_else(corrected_flag, "low", as.character(reviewer_risk_level)),
      category = recode(
        explanation_category_plot,
        disease_rescue_model = "Disease-rescue model",
        distal_or_adaptive_tissue = "Distal/adaptive tissue context",
        late_timepoint = "Late timepoint",
        formulation_designed_to_reduce_inflammation = "Low-inflammatory formulation design",
        therapeutic_cargo_specific_effect = "Therapeutic cargo/context effect",
        tissue_time_kinetic_effect = "Tissue/time kinetic effect",
        .default = NA_character_
      ),
      category = factor(category, levels = category_levels),
      support_level = factor(support_level, levels = support_levels)
    )
  if (nrow(row_tbl) != 15L || any(is.na(row_tbl$category)) || any(is.na(row_tbl$support_level))) {
    stop("Supplementary Figure S1B requires 15 fully classified context rows.", call. = FALSE)
  }

  plot_tbl <- row_tbl %>%
    count(category, support_level, name = "n", .drop = FALSE)
  total_tbl <- row_tbl %>%
    count(category, name = "total", .drop = FALSE)
  expected_totals <- c(2L, 1L, 1L, 1L, 4L, 6L)
  if (!identical(as.integer(total_tbl$total), expected_totals)) {
    stop("Supplementary Figure S1B category totals do not match 2, 1, 1, 1, 4, 6.",
         call. = FALSE)
  }
  expected_support <- matrix(
    c(0L, 2L, 0L,
      1L, 0L, 0L,
      0L, 1L, 0L,
      1L, 0L, 0L,
      1L, 3L, 0L,
      6L, 0L, 0L),
    nrow = length(category_levels), byrow = TRUE
  )
  observed_support <- matrix(
    as.integer(plot_tbl$n),
    nrow = length(category_levels), byrow = TRUE
  )
  if (!identical(observed_support, expected_support)) {
    stop("Supplementary Figure S1B support-level stacks do not match the corrected specification.",
         call. = FALSE)
  }
  segment_label_tbl <- plot_tbl %>%
    filter(category == "Therapeutic cargo/context effect", n > 0L)

  p <- ggplot(plot_tbl, aes(x = category, y = n, fill = support_level)) +
    geom_col(width = 0.7) +
    geom_text(
      data = segment_label_tbl,
      aes(label = n),
      position = position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3.1
    ) +
    geom_text(
      data = total_tbl,
      aes(x = category, y = total, label = total),
      inherit.aes = FALSE,
      vjust = -0.45,
      color = "#111827",
      fontface = "bold",
      size = 3.2
    ) +
    scale_fill_manual(
      values = support_palette,
      limits = support_levels,
      breaks = support_levels,
      labels = support_levels,
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = 0:6,
      expand = expansion(mult = c(0, 0.12))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Weak responses are explained by timing and biological context",
      subtitle = "Fifteen context-shifted contrasts are grouped by biological interpretation and support level.",
      x = "Biological interpretation category",
      y = "Number of contrasts",
      fill = "Interpretation support level"
    ) +
    theme_imrs_publication(base_size = 9.5) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.box.just = "center",
      plot.margin = margin(10, 16, 18, 10)
    )
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1B_weak_dataset_context_summary", 9, 5.5, dpi = 400,
                 source_tables = required_paths$weak_context,
                 source_code_section_or_function = "make_FigureSB_simplified")
}

make_FigureSC <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) != "Locked anchor") %>%
    group_by(dataset_id, tissue, time_h, delivery_platform_clean, manuscript_group) %>%
    summarise(mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              n_contrasts = n(), .groups = "drop") %>%
    mutate(label = short_text(dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE), 55),
           label = ordered_factor(label, mean_delta))
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_point(alpha = 0.95) +
    scale_color_manual(values = manuscript_group_palette, drop = FALSE) +
    scale_size_continuous(range = c(2, 5)) +
    labs(
      title = "Dataset-level summaries clarify context-dependent IMRS responses",
      subtitle = "Dataset-level means reduce contrast-level crowding; the full forest is retained in the appendix.",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "Dataset context",
      color = "Validation group",
      size = "Contrasts"
    ) +
    theme_imrs_publication(base_size = 9)
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1C_simplified_by_dataset", 8.2, 5.5, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSC")
}

make_FigureSD <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% c("Extended validation", "Secondary support")) %>%
    group_by(dataset_id, tissue, delivery_platform_clean, manuscript_group) %>%
    summarise(mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              n_contrasts = n(),
              time_min = min(time_h, na.rm = TRUE),
              time_max = max(time_h, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(label = short_text(dataset_context_label(dataset_id, tissue, time_min, delivery_platform_clean, compact = TRUE), 64),
           label = ordered_factor(label, mean_delta))
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_point(alpha = 0.95) +
    coord_cartesian(xlim = c(-2, 5)) +
    scale_color_manual(values = manuscript_group_palette,
                       breaks = c("Extended validation", "Secondary support"),
                       drop = FALSE) +
    scale_size_continuous(range = c(2.2, 5)) +
    guides(
      size = guide_legend(order = 1, nrow = 1),
      color = guide_legend(order = 2, nrow = 1)
    ) +
    labs(
      title = "Weak-context datasets show attenuated IMRS responses",
      subtitle = "Dataset-level means highlight late or context-shifted validation settings.",
      x = "Mean delivery-minus-control IMRS z-score",
      y = "Dataset context",
      color = "Analysis group",
      size = "Contrasts"
    ) +
    theme_imrs_publication(base_size = 8.5) +
    theme(plot.margin = margin(8, 18, 8, 10))
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1D_weak_zoom_forest", 8.8, 4.9, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSD")
}

make_FigureSE <- function() {
  plot_tbl <- weak_tbl %>%
    mutate(
      original_IMRS_delta = safe_num(original_IMRS_delta),
      context_label = dataset_context_label(dataset_id, tissue, safe_num(time_h), delivery_platform, compact = TRUE),
      explanation_label = display_text(explanation_category),
      time_label = format_time_label(time_h)
    )
  p <- ggplot(plot_tbl, aes(x = context_label, y = original_IMRS_delta, color = explanation_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563") +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_text(aes(label = time_label), vjust = -0.75, size = 2.5, show.legend = FALSE) +
    labs(
      title = "Context-shifted datasets show variable IMRS elevation",
      subtitle = "Points show delivery-minus-control IMRS shifts annotated by collection time.",
      x = "Dataset and biological context",
      y = "Mean delivery-minus-control IMRS z-score",
      color = "Explanation category"
    ) +
    theme_imrs_publication(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3A_weak_context_plot", 7.8, 4.8, dpi = 400,
                 source_tables = required_paths$weak_context,
                 source_code_section_or_function = "make_FigureSE")
}

make_FigureSF <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(dataset_id == "GSE264344") %>%
    mutate(
      tissue_label = case_when(
        str_to_lower(tissue) == "dln" ~ "dLN (draining lymph node)",
        str_to_lower(tissue) == "draining lymph nodes" ~ "dLN (draining lymph node)",
        TRUE ~ display_text(tissue)
      ),
      vector_label = case_when(
        str_detect(str_to_lower(split_id), "ad26") ~ "Ad26",
        str_detect(str_to_lower(split_id), "ad5") ~ "Ad5",
        TRUE ~ short_text(display_text(delivery_platform_clean), 18)
      ),
      point_label = ifelse(time_h == 72, "72 h", "")
    )
  p <- ggplot(plot_tbl, aes(x = time_h, y = delta_mean_imrs_z,
                            color = tissue_label, shape = vector_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_vline(xintercept = 24, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_line(aes(group = interaction(tissue_label, vector_label)), alpha = 0.55, linewidth = 0.45) +
    geom_point(size = 2.4, alpha = 0.92) +
    geom_text(aes(label = point_label), nudge_y = 0.45, size = 2.4, show.legend = FALSE) +
    facet_wrap(~ tissue_label, nrow = 1) +
    scale_x_continuous(breaks = sort(unique(plot_tbl$time_h))) +
    labs(
      title = "GSE264344 captures adenoviral-vector IMRS kinetics",
      subtitle = "dLN means draining lymph node; 72 h values are interpreted as waning kinetics.",
      x = "Time after delivery (h)",
      y = "Mean delivery-minus-control IMRS z-score",
      color = "Tissue",
      shape = "Vector"
    ) +
    theme_imrs_publication(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_imrs_plot(p, folder_path("FigureS3_context_timecourse_and_dominance_appendix"),
                 "FigureS3B_gse264344_time_course", 7.5, 4.6, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSF")
}

# -------------------------------------------------------------------------
# D. Run section
# -------------------------------------------------------------------------

write_dataset_context_companion()

panel_functions <- list(
  Figure1A = make_Figure1A,
  Figure1B = make_Figure1B,
  Figure1C = make_Figure1C,
  Figure1D = make_Figure1D,
  Figure1D_dataset_role_summary = make_dataset_role_summary_plot,
  Figure2A = make_Figure2A,
  Figure2B = make_Figure2B,
  Figure2C = make_Figure2C,
  Figure2D = make_Figure2D,
  Figure3A = make_Figure3A,
  Figure3B = make_Figure3B,
  Figure3C = make_Figure3C,
  Figure4A = make_Figure3D_simplified,
  Figure5A = make_Figure4A,
  Figure5B = make_Figure4B,
  Figure5C = make_Figure4C,
  Figure6A = make_Figure4D,
  Figure6B = make_Figure4E,
  Figure6C = make_Figure4F,
  Figure6D = make_Figure4G,
  Figure7A = make_Figure5A,
  Figure7B = make_Figure5B,
  Figure7C = make_Figure5C,
  Figure7D = make_Figure5D,
  Figure7E = make_Figure5E,
  Figure8A = make_Figure5F,
  Figure8B = make_Figure5G,
  Figure8C = make_Figure5I,
  FigureS1A = make_FigureSA,
  FigureS1B = make_FigureSB_simplified,
  FigureS1C = make_FigureSC,
  FigureS1D = make_FigureSD,
  FigureS2A = make_Figure3D,
  FigureS2B = make_FigureSB,
  FigureS3A = make_FigureSE,
  FigureS3B = make_FigureSF,
  FigureS3C = make_Figure5H,
  FigureS3D = make_Figure2E,
  FigureS3E = make_Figure3E,
  FigureS3F = make_Figure3F
)

for (panel_name in names(panel_functions)) {
  fn <- panel_functions[[panel_name]]
  tryCatch(
    fn(),
    error = function(e) {
      skipped_panels <<- bind_rows(skipped_panels, tibble(panel_id = panel_name, reason = conditionMessage(e)))
      write_log()
      stop("Figure generation failed in ", panel_name, ": ", conditionMessage(e), call. = FALSE)
    }
  )
}

manifest <- bind_rows(manifest_rows)
validate_manifest(manifest)

readr::write_tsv(manifest, manifest_path, na = "NA")

log_msg("Wrote manifest: ", manifest_path)
log_msg("Completed IMRS figure generation.")
write_log()

total_png <- sum(file.exists(manifest$output_png))
total_pdf <- sum(!is.na(manifest$output_pdf) & file.exists(manifest$output_pdf))
total_svg <- sum(!is.na(manifest$output_svg) & file.exists(manifest$output_svg))
combined_figures_generated <- 0L

cat("\nIMRS final figure regeneration summary\n")
cat("--------------------------------------\n")
cat("Standalone panels generated: ", nrow(manifest), "\n", sep = "")
cat("Combined figures generated: ", combined_figures_generated, "\n", sep = "")
cat("Total PNG files generated: ", total_png, "\n", sep = "")
cat("Total PDF files generated: ", total_pdf, "\n", sep = "")
cat("Total SVG files generated: ", total_svg, "\n", sep = "")
cat("Output root path: ", normalizePath(output_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Skipped panels: ", ifelse(nrow(skipped_panels) == 0, "none", paste(skipped_panels$panel_id, collapse = ", ")), "\n", sep = "")
cat("Warnings: ", ifelse(length(warnings_seen) == 0, "none", paste(warnings_seen, collapse = " | ")), "\n", sep = "")
cat("Manifest path: ", normalizePath(manifest_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Log path: ", normalizePath(log_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
