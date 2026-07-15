#!/usr/bin/env Rscript

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "03_leave_one_gene_out"

select_leave_one_genes <- function(weights_df) {
  required_symbols <- c("Acod1", "Cxcl11", "Ccl2", "Cxcl10", "Ifit1", "Ifit3", "Mx2", "Saa3")
  top <- weights_df %>%
    arrange(desc(abs(beta_meta))) %>%
    slice_head(n = DEFAULT_TOP_GENE_N) %>%
    mutate(selection_reason = "top_abs_beta_meta")
  required <- weights_df %>%
    filter(tolower(gene_symbol) %in% tolower(required_symbols)) %>%
    mutate(selection_reason = "specified_gene_symbol")
  bind_rows(top, required) %>%
    arrange(desc(abs(beta_meta))) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(removal_rank = row_number())
}

score_removed_gene <- function(comp, remove_gene_id) {
  keep <- setdiff(names(comp$wvec), remove_gene_id)
  score_from_z(comp$z[keep, , drop = FALSE], comp$wvec[keep], comp$design,
               comp$coverage, DEFAULT_SCORE_SD_FLOOR)
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  weights <- load_frozen_weights()
  remove_genes <- select_leave_one_genes(weights)
  if (nrow(remove_genes) == 0) stop("No genes selected for leave-one-gene-out.")

  eval_tbl <- load_step09_eval(required = TRUE) %>%
    mutate(pass = as.logical(pass)) %>%
    filter(pass %in% TRUE)
  if (nrow(eval_tbl) == 0) stop("No passing Step 09 contrasts available.")

  datasets <- sort(unique(eval_tbl$gse_id))
  log_msg("INFO", "Scoring ", length(datasets), " datasets for leave-one-gene-out.")

  all_rows <- map_dfr(datasets, function(ds) {
    rows <- eval_tbl %>% filter(gse_id == ds)
    comp <- tryCatch(
      score_dataset_components(ds, weights, enforce_min_coverage = TRUE),
      error = function(e) {
        log_msg("WARN", "Skipping leave-one-gene-out dataset ", ds, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(comp)) return(tibble())

    original_eval <- evaluate_scores_for_step09_rows(comp$score_tbl, rows) %>%
      select(split_path, original_delta_mean_imrs_z = delta_mean_imrs_z,
             original_auc_imrs_z = auc_imrs_z)

    map_dfr(seq_len(nrow(remove_genes)), function(i) {
      gene <- remove_genes[i, ]
      score_tbl <- if (gene$gene_id %in% names(comp$wvec)) {
        score_removed_gene(comp, gene$gene_id)
      } else {
        comp$score_tbl
      }
      loo_eval <- evaluate_scores_for_step09_rows(score_tbl, rows) %>%
        left_join(original_eval, by = "split_path") %>%
        mutate(
          removed_gene_id = gene$gene_id,
          removed_gene_symbol = gene$gene_symbol,
          removed_gene_beta_meta = gene$beta_meta,
          removal_rank = gene$removal_rank,
          removal_selection_reason = gene$selection_reason,
          gene_present_in_dataset = gene$gene_id %in% names(comp$wvec),
          leave_one_gene_out_delta_mean_imrs_z = delta_mean_imrs_z,
          absolute_change_delta = abs(leave_one_gene_out_delta_mean_imrs_z -
                                        original_delta_mean_imrs_z),
          percent_change_delta = ifelse(
            is.finite(original_delta_mean_imrs_z) & abs(original_delta_mean_imrs_z) > 0,
            100 * (leave_one_gene_out_delta_mean_imrs_z - original_delta_mean_imrs_z) /
              abs(original_delta_mean_imrs_z),
            NA_real_
          ),
          absolute_percent_change_delta = abs(percent_change_delta),
          direction_preserved = (original_delta_mean_imrs_z > 0) ==
            (leave_one_gene_out_delta_mean_imrs_z > 0),
          audit_role = assign_audit_role(gse_id, time_h),
          time_group = assign_time_group(time_h)
        ) %>%
        select(
          removed_gene_id, removed_gene_symbol, removal_rank,
          removal_selection_reason, removed_gene_beta_meta, gene_present_in_dataset,
          gse_id, split_id, split_path, audit_role, time_group, tissue, time_h,
          n_controls, n_delivery,
          original_delta_mean_imrs_z, leave_one_gene_out_delta_mean_imrs_z,
          absolute_change_delta, percent_change_delta, absolute_percent_change_delta,
          direction_preserved, original_auc_imrs_z, auc_imrs_z
        )
      loo_eval
    })
  })

  out_path <- file.path(AUDIT_RESULTS_DIR, "leave_one_gene_out_summary.tsv")
  write_tsv_safe(all_rows, out_path)

  if (nrow(all_rows) == 0) {
    p_corr <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No leave-one-gene-out results") +
      theme_void()
    p_effect <- p_corr
  } else {
    p_corr <- ggplot(all_rows, aes(x = original_delta_mean_imrs_z,
                                   y = leave_one_gene_out_delta_mean_imrs_z)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
      geom_point(aes(color = audit_role), alpha = 0.45, size = 1.4) +
      scale_color_manual(values = c(anchor = "#4E79A7", calibration = "#59A14F",
                                    external = "#E15759"), drop = FALSE) +
      labs(
        title = "Leave-One-Gene-Out Delta Correlation",
        x = "Original delta IMRSz",
        y = "Leave-one-gene-out delta IMRSz",
        color = "Audit role"
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))

    effect_tbl <- all_rows %>%
      group_by(removed_gene_id, removed_gene_symbol, removal_rank) %>%
      summarise(
        max_absolute_change_delta = max(absolute_change_delta, na.rm = TRUE),
        median_absolute_percent_change_delta = safe_median(absolute_percent_change_delta),
        n_direction_not_preserved = sum(!direction_preserved, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(max_absolute_change_delta)) %>%
      slice_head(n = 20) %>%
      mutate(label = ifelse(!is.na(removed_gene_symbol) & nzchar(removed_gene_symbol),
                            removed_gene_symbol, removed_gene_id),
             label = factor(label, levels = rev(unique(label))))

    p_effect <- ggplot(effect_tbl, aes(x = max_absolute_change_delta, y = label,
                                       fill = n_direction_not_preserved > 0)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#E15759"),
                        labels = c(`FALSE` = "Direction preserved", `TRUE` = "Any direction change")) +
      labs(
        title = "Largest Delta Change After Removing One Gene",
        x = "Maximum absolute delta change",
        y = NULL,
        fill = NULL
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  save_plot(p_corr, file.path(AUDIT_PLOTS_DIR,
                              "leave_one_gene_out_delta_correlation.png"),
            width = 7, height = 5)
  save_plot(p_effect, file.path(AUDIT_PLOTS_DIR, "top_gene_removal_effect.png"),
            width = 7.5, height = 5.5)

  log_msg("INFO", "Wrote leave-one-gene-out rows=", nrow(all_rows), " to ", out_path)
  invisible(all_rows)
}

main()
