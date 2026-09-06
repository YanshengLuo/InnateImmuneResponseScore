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

this_script <- imrs_detect_script_path_local("05_threshold_sensitivity.R")
if (is.na(this_script)) stop("Could not identify this script path. Use RStudio Source/Run or Rscript.", call. = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "05_threshold_sensitivity"

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  eval_tbl <- load_step09_eval(required = TRUE) %>%
    mutate(pass = as.logical(pass),
           audit_role = assign_audit_role(gse_id, time_h),
           time_group = assign_time_group(time_h)) %>%
    filter(pass %in% TRUE)
  if (nrow(eval_tbl) == 0) stop("No passing Step 09 contrasts available.")

  frozen_genes <- load_frozen_weights() %>% pull(gene_id) %>% unique()
  grid <- expand_grid(
    min_abs_log2FC = c(0.25, 0.50, 0.75),
    min_up_anchor_count = c(2L, 3L),
    fdr_support = c(0.05, 0.10, 0.20)
  ) %>%
    mutate(grid_id = paste0("lfc", min_abs_log2FC, "_k", min_up_anchor_count,
                            "_fdr", fdr_support))

  external_rows <- eval_tbl %>% filter(audit_role != "anchor")
  anchor_rows <- eval_tbl %>% filter(audit_role == "anchor", gse_id %in% STRICT_ANCHOR_IDS)
  needed_datasets <- sort(unique(c(external_rows$gse_id, anchor_rows$gse_id)))

  prep_cache <- new.env(parent = emptyenv())
  get_prep <- function(ds) {
    if (!exists(ds, envir = prep_cache, inherits = FALSE)) {
      log_msg("INFO", "Preparing expression cache for ", ds)
      assign(ds, prepare_dataset_expression(ds), envir = prep_cache)
    }
    get(ds, envir = prep_cache, inherits = FALSE)
  }

  detail_rows <- list()
  gene_set_rows <- list()
  summary_tbl <- map_dfr(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    log_msg("INFO", "Threshold grid ", i, "/", nrow(grid), ": ", g$grid_id)

    full_weights <- tryCatch(
      build_weights_from_de(
        training_anchor_ids = STRICT_ANCHOR_IDS,
        min_abs_log2fc = g$min_abs_log2FC,
        fdr_support = g$fdr_support,
        min_up_anchor_count = g$min_up_anchor_count,
        require_up = TRUE
      ),
      error = function(e) {
        log_msg("WARN", "Full threshold weight build failed for ", g$grid_id,
                ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(full_weights) || nrow(full_weights) == 0) {
      return(tibble(
        grid_id = g$grid_id,
        min_abs_log2FC = g$min_abs_log2FC,
        min_up_anchor_count = g$min_up_anchor_count,
        fdr_support = g$fdr_support,
        n_genes_selected = 0L,
        overlap_with_frozen_genes = 0L,
        overlap_fraction_with_frozen = NA_real_,
        external_n_contrasts = 0L,
        external_mean_delta_imrs_z = NA_real_,
        external_median_delta_imrs_z = NA_real_,
        external_proportion_positive_delta = NA_real_,
        external_mean_auc_secondary = NA_real_,
        external_median_auc_secondary = NA_real_,
        loao_n_contrasts = 0L,
        loao_mean_delta_imrs_z = NA_real_,
        loao_proportion_positive_delta = NA_real_,
        status = "failed_weight_build"
      ))
    }

    gene_set_rows[[g$grid_id]] <<- full_weights %>%
      transmute(
        grid_id = g$grid_id,
        min_abs_log2FC = g$min_abs_log2FC,
        min_up_anchor_count = g$min_up_anchor_count,
        fdr_support = g$fdr_support,
        gene_id = gene_id,
        beta_meta = beta_meta,
        se_meta = se_meta,
        I2 = I2,
        support_datasets = support_datasets,
        datasets_supporting = datasets_supporting,
        in_frozen_gene_set = gene_id %in% frozen_genes
      )

    external_eval <- map_dfr(sort(unique(external_rows$gse_id)), function(ds) {
      rows <- external_rows %>% filter(gse_id == ds)
      comp <- tryCatch(
        score_prepared_dataset(get_prep(ds), full_weights, enforce_min_coverage = FALSE),
        error = function(e) {
          log_msg("WARN", "Threshold scoring failed for ", ds, " / ", g$grid_id,
                  ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(comp)) return(tibble())
      evaluate_scores_for_step09_rows(comp$score_tbl, rows) %>%
        mutate(
          grid_id = g$grid_id,
          min_abs_log2FC = g$min_abs_log2FC,
          min_up_anchor_count = g$min_up_anchor_count,
          fdr_support = g$fdr_support,
          audit_role = assign_audit_role(gse_id, time_h),
          time_group = assign_time_group(time_h),
          sensitivity_scope = "external_full_3_anchor_weights"
        )
    })

    loao_eval <- map_dfr(STRICT_ANCHOR_IDS, function(holdout) {
      training <- setdiff(STRICT_ANCHOR_IDS, holdout)
      rows <- anchor_rows %>% filter(gse_id == holdout)
      if (nrow(rows) == 0) return(tibble())
      w2 <- tryCatch(
        build_weights_from_de(
          training_anchor_ids = training,
          min_abs_log2fc = g$min_abs_log2FC,
          fdr_support = g$fdr_support,
          min_up_anchor_count = length(training),
          require_up = TRUE
        ),
        error = function(e) {
          log_msg("WARN", "LOAO threshold weight build failed for holdout ",
                  holdout, " / ", g$grid_id, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(w2) || nrow(w2) == 0) return(tibble())
      comp <- tryCatch(
        score_prepared_dataset(get_prep(holdout), w2, enforce_min_coverage = FALSE),
        error = function(e) {
          log_msg("WARN", "LOAO threshold scoring failed for holdout ", holdout,
                  " / ", g$grid_id, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(comp)) return(tibble())
      evaluate_scores_for_step09_rows(comp$score_tbl, rows) %>%
        mutate(
          grid_id = g$grid_id,
          min_abs_log2FC = g$min_abs_log2FC,
          min_up_anchor_count = g$min_up_anchor_count,
          fdr_support = g$fdr_support,
          holdout_anchor = holdout,
          training_anchors = paste(training, collapse = ";"),
          n_training_genes = nrow(w2),
          audit_role = "anchor",
          time_group = assign_time_group(time_h),
          sensitivity_scope = "loao_2_anchor_weights"
        )
    })

    detail_rows[[g$grid_id]] <<- bind_rows(external_eval, loao_eval)

    tibble(
      grid_id = g$grid_id,
      min_abs_log2FC = g$min_abs_log2FC,
      min_up_anchor_count = g$min_up_anchor_count,
      fdr_support = g$fdr_support,
      n_genes_selected = nrow(full_weights),
      overlap_with_frozen_genes = sum(full_weights$gene_id %in% frozen_genes),
      overlap_fraction_with_frozen = sum(full_weights$gene_id %in% frozen_genes) /
        length(frozen_genes),
      external_n_contrasts = nrow(external_eval),
      external_mean_delta_imrs_z = safe_mean(external_eval$delta_mean_imrs_z),
      external_median_delta_imrs_z = safe_median(external_eval$delta_mean_imrs_z),
      external_proportion_positive_delta = safe_prop(external_eval$delta_mean_imrs_z > 0),
      external_mean_auc_secondary = safe_mean(external_eval$auc_imrs_z),
      external_median_auc_secondary = safe_median(external_eval$auc_imrs_z),
      loao_n_contrasts = nrow(loao_eval),
      loao_mean_delta_imrs_z = safe_mean(loao_eval$delta_mean_imrs_z),
      loao_proportion_positive_delta = safe_prop(loao_eval$delta_mean_imrs_z > 0),
      status = "ok"
    )
  })

  gene_sets <- bind_rows(gene_set_rows)
  detail_tbl <- bind_rows(detail_rows)

  out_summary <- file.path(AUDIT_RESULTS_DIR, "threshold_sensitivity_summary.tsv")
  out_genes <- file.path(AUDIT_RESULTS_DIR, "threshold_gene_sets.tsv")
  out_detail <- file.path(AUDIT_RESULTS_DIR, "threshold_sensitivity_contrast_deltas.tsv")
  write_tsv_safe(summary_tbl, out_summary)
  write_tsv_safe(gene_sets, out_genes)
  write_tsv_safe(detail_tbl, out_detail)

  if (nrow(summary_tbl) == 0 || all(!is.finite(summary_tbl$external_mean_delta_imrs_z))) {
    p_heat <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No threshold sensitivity results") +
      theme_void()
  } else {
    p_heat <- ggplot(summary_tbl,
                     aes(x = factor(min_abs_log2FC), y = factor(fdr_support),
                         fill = external_mean_delta_imrs_z)) +
      geom_tile(color = "white", linewidth = 0.6) +
      geom_text(aes(label = paste0("n=", n_genes_selected, "\n",
                                   sprintf("%.2f", external_mean_delta_imrs_z))),
                size = 3) +
      facet_wrap(~ min_up_anchor_count, labeller = label_both) +
      scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                           midpoint = 0, na.value = "grey85") +
      labs(
        title = "Threshold Sensitivity: External Mean Delta",
        x = "Minimum log2FC",
        y = "FDR support",
        fill = "Mean delta"
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  dist_tbl <- detail_tbl %>%
    filter(sensitivity_scope == "external_full_3_anchor_weights")
  if (nrow(dist_tbl) == 0) {
    p_dist <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No external threshold deltas") +
      theme_void()
  } else {
    p_dist <- ggplot(dist_tbl,
                     aes(x = factor(min_abs_log2FC), y = delta_mean_imrs_z,
                         fill = factor(fdr_support))) +
      geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
      geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
      facet_wrap(~ min_up_anchor_count, labeller = label_both) +
      scale_fill_manual(values = c("0.05" = "#4E79A7", "0.1" = "#59A14F",
                                   "0.2" = "#F28E2B"), drop = FALSE) +
      labs(
        title = "Threshold Sensitivity: External Delta Distribution",
        x = "Minimum log2FC",
        y = "Delta IMRSz",
        fill = "FDR"
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  save_plot(p_heat, file.path(AUDIT_PLOTS_DIR, "threshold_sensitivity_heatmap.png"),
            width = 8.5, height = 5.5)
  save_plot(p_dist, file.path(AUDIT_PLOTS_DIR,
                              "threshold_sensitivity_delta_distribution.png"),
            width = 8.5, height = 5.5)

  log_msg("INFO", "Wrote threshold sensitivity rows=", nrow(summary_tbl),
          " to ", out_summary)
  invisible(summary_tbl)
}

main()
