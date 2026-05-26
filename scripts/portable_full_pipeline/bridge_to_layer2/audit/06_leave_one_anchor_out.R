#!/usr/bin/env Rscript

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "06_leave_one_anchor_out"

collect_split_samples <- function(score_tbl, row, holdout_anchor, training_anchors) {
  split <- read_split_table(row$split_path)
  split %>%
    select(any_of(c("sample_id", "condition_simple", "tissue", "time_h",
                    "contrast_label", "control_label"))) %>%
    rename(split_condition_simple = condition_simple) %>%
    inner_join(score_tbl %>% select(sample_id, imrs_raw, imrs_z, coverage),
               by = "sample_id") %>%
    mutate(
      condition_simple = split_condition_simple,
      holdout_anchor = holdout_anchor,
      training_anchors = paste(training_anchors, collapse = ";"),
      split_id = row$split_id,
      gse_id = row$gse_id,
      tissue = if ("tissue" %in% names(.)) tissue else row$tissue,
      time_h = suppressWarnings(as.numeric(if ("time_h" %in% names(.)) time_h else row$time_h))
    )
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  eval_tbl <- load_step09_eval(required = TRUE) %>%
    mutate(pass = as.logical(pass),
           audit_role = assign_audit_role(gse_id, time_h)) %>%
    filter(pass %in% TRUE, audit_role == "anchor", gse_id %in% STRICT_ANCHOR_IDS)
  if (nrow(eval_tbl) == 0) stop("No strict-anchor Step 09 contrasts available.")

  prep_cache <- new.env(parent = emptyenv())
  get_prep <- function(ds) {
    if (!exists(ds, envir = prep_cache, inherits = FALSE)) {
      assign(ds, prepare_dataset_expression(ds), envir = prep_cache)
    }
    get(ds, envir = prep_cache, inherits = FALSE)
  }

  sample_list <- list()
  detail_tbl <- map_dfr(STRICT_ANCHOR_IDS, function(holdout) {
    training <- setdiff(STRICT_ANCHOR_IDS, holdout)
    rows <- eval_tbl %>% filter(gse_id == holdout)
    if (nrow(rows) == 0) {
      log_msg("WARN", "No held-out split rows for strict anchor ", holdout)
      return(tibble())
    }
    weights <- tryCatch(
      build_weights_from_de(
        training_anchor_ids = training,
        min_abs_log2fc = 0.50,
        fdr_support = 0.05,
        min_up_anchor_count = length(training),
        require_up = TRUE
      ),
      error = function(e) {
        log_msg("WARN", "LOAO weight build failed for holdout ", holdout, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(weights) || nrow(weights) == 0) return(tibble())
    comp <- tryCatch(
      score_prepared_dataset(get_prep(holdout), weights, enforce_min_coverage = FALSE),
      error = function(e) {
        log_msg("WARN", "LOAO scoring failed for holdout ", holdout, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(comp)) return(tibble())
    sample_rows <- map_dfr(seq_len(nrow(rows)), function(i) {
      tryCatch(
        collect_split_samples(comp$score_tbl, rows[i, ], holdout, training),
        error = function(e) tibble()
      )
    })
    if (nrow(sample_rows) > 0) sample_list[[holdout]] <<- sample_rows
    evaluate_scores_for_step09_rows(comp$score_tbl, rows) %>%
      mutate(
        holdout_anchor = holdout,
        training_anchors = paste(training, collapse = ";"),
        n_training_genes = nrow(weights),
        n_training_gene_overlap_heldout = comp$n_overlap,
        heldout_coverage = comp$coverage,
        direction_positive = delta_mean_imrs_z > 0,
        auc_secondary = auc_imrs_z
      )
  })

  sample_tbl <- bind_rows(sample_list)

  summary_tbl <- detail_tbl %>%
    group_by(holdout_anchor, training_anchors) %>%
    summarise(
      n_training_genes = first(n_training_genes),
      n_training_gene_overlap_heldout = first(n_training_gene_overlap_heldout),
      heldout_coverage = first(heldout_coverage),
      n_contrasts = n(),
      mean_delta_imrs_z = safe_mean(delta_mean_imrs_z),
      median_delta_imrs_z = safe_median(delta_mean_imrs_z),
      proportion_positive_delta = safe_prop(delta_mean_imrs_z > 0),
      all_directions_positive = all(delta_mean_imrs_z > 0, na.rm = TRUE),
      mean_auc_secondary = safe_mean(auc_secondary),
      median_auc_secondary = safe_median(auc_secondary),
      .groups = "drop"
    ) %>%
    arrange(match(holdout_anchor, STRICT_ANCHOR_IDS))

  out_summary <- file.path(AUDIT_RESULTS_DIR, "leave_one_anchor_out_summary.tsv")
  out_detail <- file.path(AUDIT_RESULTS_DIR, "leave_one_anchor_out_contrast_details.tsv")
  write_tsv_safe(summary_tbl, out_summary)
  write_tsv_safe(detail_tbl, out_detail)

  if (nrow(detail_tbl) == 0) {
    p_delta <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No leave-one-anchor-out results") +
      theme_void()
  } else {
    p_delta <- ggplot(detail_tbl, aes(x = holdout_anchor, y = delta_mean_imrs_z,
                                      color = holdout_anchor)) +
      geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
      geom_jitter(width = 0.16, height = 0, size = 2, alpha = 0.75) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.45,
                   color = "black", linewidth = 0.35) +
      scale_color_manual(values = c(GSE264344 = "#4E79A7", GSE39129 = "#59A14F",
                                    GSE167521 = "#E15759"), drop = FALSE) +
      labs(
        title = "Leave-One-Anchor-Out Delta",
        x = "Held-out strict anchor",
        y = "Delta IMRSz (DELIVERY - CONTROL)",
        color = "Held-out anchor"
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  if (nrow(sample_tbl) == 0) {
    p_samples <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No held-out sample scores") +
      theme_void()
  } else {
    p_samples <- ggplot(sample_tbl, aes(x = condition_simple, y = imrs_z,
                                        fill = condition_simple)) +
      geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
      geom_violin(alpha = 0.55, trim = FALSE, width = 0.8) +
      geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.08, size = 0.8, alpha = 0.45) +
      facet_wrap(~ holdout_anchor, scales = "free_y") +
      scale_fill_manual(values = c(CONTROL = "#9E9E9E", DELIVERY = "#4E79A7")) +
      labs(
        title = "Held-Out Anchor Sample Score Distributions",
        x = NULL,
        y = "IMRSz",
        fill = NULL
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  save_plot(p_delta, file.path(AUDIT_PLOTS_DIR, "leave_one_anchor_out_delta.png"),
            width = 7.2, height = 5)
  save_plot(p_samples, file.path(AUDIT_PLOTS_DIR,
                                 "leave_one_anchor_out_sample_distributions.png"),
            width = 8.5, height = 5.5)

  log_msg("INFO", "Wrote leave-one-anchor-out summary rows=", nrow(summary_tbl),
          " to ", out_summary)
  invisible(summary_tbl)
}

main()
