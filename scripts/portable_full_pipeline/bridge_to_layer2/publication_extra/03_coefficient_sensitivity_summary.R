#!/usr/bin/env Rscript

# Summarize existing coefficient/sensitivity outputs without rerunning them.
# Missing files are explicitly inventoried rather than inferred.

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_publication_extra_config.R"))

existing_sensitivity_files <- tribble(
  ~analysis, ~path, ~description,
  "leave_one_gene_out",
  file.path(AUDIT_RESULTS_DIR, "leave_one_gene_out_summary.tsv"),
  "Existing single-gene coefficient removal sensitivity.",
  "leave_one_anchor_out",
  file.path(AUDIT_RESULTS_DIR, "leave_one_anchor_out_summary.tsv"),
  "Existing anchor holdout robustness check.",
  "threshold_sensitivity",
  file.path(AUDIT_RESULTS_DIR, "threshold_sensitivity_summary.tsv"),
  "Existing gene-selection threshold sensitivity grid.",
  "threshold_sensitivity_contrasts",
  file.path(AUDIT_RESULTS_DIR, "threshold_sensitivity_contrast_deltas.tsv"),
  "Existing per-contrast deltas from threshold sensitivity grid.",
  "gene_dominance",
  file.path(AUDIT_RESULTS_DIR, "gene_dominance_summary.tsv"),
  "Existing top-contributor dominance audit."
)

count_rows <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    lines <- readLines(con, n = 10000, warn = FALSE)
    if (length(lines) == 0) break
    n <- n + length(lines)
  }
  max(n - 1L, 0L)
}

add_metric <- function(metric, value, interpretation, source_file) {
  tibble(
    metric = metric,
    value = as.character(value),
    interpretation = interpretation,
    source_file = source_file
  )
}

main <- function() {
  log_extra("Starting coefficient sensitivity summary from existing outputs.")

  inventory <- existing_sensitivity_files %>%
    mutate(
      exists = file.exists(path),
      row_count = map_int(path, count_rows),
      size_bytes = ifelse(exists, file.info(path)$size, NA_real_),
      modified_time = ifelse(exists, as.character(file.info(path)$mtime),
                             NA_character_)
    )

  metrics <- list()

  logo_path <- inventory$path[inventory$analysis == "leave_one_gene_out"]
  if (file.exists(logo_path)) {
    logo <- readr::read_tsv(logo_path, show_col_types = FALSE, progress = FALSE)
    metrics <- append(metrics, list(
      add_metric(
        "leave_one_gene_out_rows", nrow(logo),
        "Rows are gene-removal by split-contrast records from the existing coefficient sensitivity run.",
        logo_path
      ),
      add_metric(
        "leave_one_gene_out_removed_genes",
        dplyr::n_distinct(logo$removed_gene_id),
        "Number of unique removed genes represented in the existing sensitivity output.",
        logo_path
      ),
      add_metric(
        "leave_one_gene_out_direction_preserved_fraction",
        round(safe_prop(as.logical(logo$direction_preserved)), 4),
        "Fraction of gene-removal contrast records where Delta IMRSz direction stayed the same.",
        logo_path
      ),
      add_metric(
        "leave_one_gene_out_median_absolute_percent_change_delta",
        round(safe_median(logo$absolute_percent_change_delta), 4),
        "Median absolute percent change in Delta IMRSz after single-gene removal.",
        logo_path
      ),
      add_metric(
        "leave_one_gene_out_max_absolute_delta_change",
        round(max(suppressWarnings(as.numeric(logo$absolute_change_delta)),
                  na.rm = TRUE), 4),
        "Largest observed absolute Delta IMRSz change after removing one gene.",
        logo_path
      )
    ))
  }

  loao_path <- inventory$path[inventory$analysis == "leave_one_anchor_out"]
  if (file.exists(loao_path)) {
    loao <- readr::read_tsv(loao_path, show_col_types = FALSE, progress = FALSE)
    metrics <- append(metrics, list(
      add_metric(
        "leave_one_anchor_out_holdouts", nrow(loao),
        "Number of strict-anchor holdout summaries in the existing robustness output.",
        loao_path
      ),
      add_metric(
        "leave_one_anchor_out_all_positive_holdouts",
        all(as.logical(loao$all_directions_positive), na.rm = TRUE),
        "Whether every held-out strict-anchor summary retained positive split-level directionality.",
        loao_path
      ),
      add_metric(
        "leave_one_anchor_out_min_mean_delta",
        round(min(suppressWarnings(as.numeric(loao$mean_delta_imrs_z)),
                  na.rm = TRUE), 4),
        "Smallest mean Delta IMRSz across strict-anchor holdout summaries.",
        loao_path
      )
    ))
  }

  threshold_path <- inventory$path[inventory$analysis == "threshold_sensitivity"]
  if (file.exists(threshold_path)) {
    threshold <- readr::read_tsv(threshold_path, show_col_types = FALSE,
                                 progress = FALSE)
    metrics <- append(metrics, list(
      add_metric(
        "threshold_sensitivity_grid_count", nrow(threshold),
        "Number of coefficient/gene-selection threshold settings already evaluated.",
        threshold_path
      ),
      add_metric(
        "threshold_sensitivity_external_positive_delta_min",
        round(min(suppressWarnings(as.numeric(
          threshold$external_proportion_positive_delta)), na.rm = TRUE), 4),
        "Minimum positive-delta fraction across external contrasts in the threshold grid.",
        threshold_path
      ),
      add_metric(
        "threshold_sensitivity_external_mean_delta_range",
        paste(
          round(range(suppressWarnings(as.numeric(
            threshold$external_mean_delta_imrs_z)), na.rm = TRUE), 4),
          collapse = " to "
        ),
        "Range of external mean Delta IMRSz across threshold settings.",
        threshold_path
      ),
      add_metric(
        "threshold_sensitivity_status_values",
        paste(sort(unique(threshold$status)), collapse = ";"),
        "Run-status values reported by the existing threshold sensitivity output.",
        threshold_path
      )
    ))
  }

  dominance_path <- inventory$path[inventory$analysis == "gene_dominance"]
  if (file.exists(dominance_path)) {
    dominance <- readr::read_tsv(dominance_path, show_col_types = FALSE,
                                 progress = FALSE)
    metrics <- append(metrics, list(
      add_metric(
        "gene_dominance_rows", nrow(dominance),
        "Number of split contrasts represented in the existing top-contributor dominance audit.",
        dominance_path
      ),
      add_metric(
        "gene_dominance_median_max_contribution_fraction",
        round(safe_median(dominance$median_max_contribution_fraction), 4),
        "Median, across contrasts, of the per-sample maximum contribution fraction.",
        dominance_path
      ),
      add_metric(
        "gene_dominance_max_max_contribution_fraction",
        round(max(suppressWarnings(as.numeric(
          dominance$max_max_contribution_fraction)), na.rm = TRUE), 4),
        "Largest observed maximum contribution fraction in the existing audit.",
        dominance_path
      )
    ))
  }

  summary_tbl <- bind_rows(metrics)
  if (nrow(summary_tbl) == 0) {
    summary_tbl <- tibble(
      metric = "no_existing_sensitivity_outputs_found",
      value = "NA",
      interpretation = paste(
        "No expected coefficient/sensitivity outputs were present.",
        "The report should flag this instead of inventing robustness evidence."
      ),
      source_file = NA_character_
    )
  }

  write_tsv_safe(inventory, file.path(EXTRA_RESULTS_DIR,
                                      "coefficient_sensitivity_inventory.tsv"))
  write_tsv_safe(summary_tbl, file.path(EXTRA_RESULTS_DIR,
                                        "coefficient_sensitivity_summary.tsv"))

  plot_tbl <- summary_tbl %>%
    filter(metric %in% c(
      "leave_one_gene_out_direction_preserved_fraction",
      "threshold_sensitivity_external_positive_delta_min",
      "gene_dominance_median_max_contribution_fraction",
      "gene_dominance_max_max_contribution_fraction"
    )) %>%
    mutate(value_num = suppressWarnings(as.numeric(value)))

  if (nrow(plot_tbl) > 0) {
    p <- ggplot(plot_tbl, aes(x = reorder(metric, value_num), y = value_num,
                              fill = metric)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      coord_flip() +
      labs(
        title = "Existing Coefficient/Sensitivity Audit Summary",
        subtitle = "Values are summarized from existing outputs; no sensitivity rerun was performed.",
        x = NULL,
        y = "Metric value"
      ) +
      theme_extra()
  } else {
    p <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = "No numeric coefficient/sensitivity outputs found") +
      theme_void()
  }
  save_plot_pair(p, "coefficient_sensitivity_existing_summary",
                 width = 7.4, height = 4.8)

  log_extra("Completed coefficient sensitivity summary. Metrics=", nrow(summary_tbl),
            ".")
  invisible(summary_tbl)
}

main()
