#!/usr/bin/env Rscript

# Baseline signature benchmarking against frozen IMRS split contrasts.
# The baseline signatures are deliberately simple and transparent; they are
# unweighted mean control-standardized gene z-scores.

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_publication_extra_config.R"))

make_signature_definitions <- function() {
  tribble(
    ~signature_id, ~signature_label, ~gene_symbol,
    "ISG", "Interferon-stimulated genes", "Ifit1",
    "ISG", "Interferon-stimulated genes", "Ifit2",
    "ISG", "Interferon-stimulated genes", "Ifit3",
    "ISG", "Interferon-stimulated genes", "Isg15",
    "ISG", "Interferon-stimulated genes", "Mx1",
    "ISG", "Interferon-stimulated genes", "Mx2",
    "ISG", "Interferon-stimulated genes", "Oas1a",
    "ISG", "Interferon-stimulated genes", "Oas1b",
    "ISG", "Interferon-stimulated genes", "Oas2",
    "ISG", "Interferon-stimulated genes", "Oas3",
    "ISG", "Interferon-stimulated genes", "Oasl1",
    "ISG", "Interferon-stimulated genes", "Oasl2",
    "ISG", "Interferon-stimulated genes", "Irf7",
    "ISG", "Interferon-stimulated genes", "Rsad2",
    "ISG", "Interferon-stimulated genes", "Usp18",
    "ISG", "Interferon-stimulated genes", "Ifih1",
    "ISG", "Interferon-stimulated genes", "Ddx58",
    "ISG", "Interferon-stimulated genes", "Stat1",
    "ISG", "Interferon-stimulated genes", "Ifi44",
    "ISG", "Interferon-stimulated genes", "Ifi44l",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ccl2",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ccl3",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ccl4",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ccl5",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ccl7",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Cxcl1",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Cxcl2",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Cxcl9",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Cxcl10",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Cxcl11",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Il1b",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Il6",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Tnf",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Nfkbia",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Ptgs2",
    "CHEMOKINE_INFLAMMATORY", "Chemokine/inflammatory genes", "Nos2",
    "GENERIC_INNATE", "Generic innate immune genes", "Tlr2",
    "GENERIC_INNATE", "Generic innate immune genes", "Tlr3",
    "GENERIC_INNATE", "Generic innate immune genes", "Tlr4",
    "GENERIC_INNATE", "Generic innate immune genes", "Tlr7",
    "GENERIC_INNATE", "Generic innate immune genes", "Tlr9",
    "GENERIC_INNATE", "Generic innate immune genes", "Myd88",
    "GENERIC_INNATE", "Generic innate immune genes", "Ticam1",
    "GENERIC_INNATE", "Generic innate immune genes", "Ddx58",
    "GENERIC_INNATE", "Generic innate immune genes", "Ifih1",
    "GENERIC_INNATE", "Generic innate immune genes", "Mavs",
    "GENERIC_INNATE", "Generic innate immune genes", "Irf3",
    "GENERIC_INNATE", "Generic innate immune genes", "Irf7",
    "GENERIC_INNATE", "Generic innate immune genes", "Nfkb1",
    "GENERIC_INNATE", "Generic innate immune genes", "Nfkbia",
    "GENERIC_INNATE", "Generic innate immune genes", "Stat1",
    "GENERIC_INNATE", "Generic innate immune genes", "Nlrp3",
    "GENERIC_INNATE", "Generic innate immune genes", "Casp1",
    "GENERIC_INNATE", "Generic innate immune genes", "Cd14",
    "GENERIC_INNATE", "Generic innate immune genes", "Ly6e"
  ) %>%
    mutate(gene_symbol_lower = tolower(gene_symbol))
}

score_signatures_for_dataset <- function(dataset_id, signature_defs,
                                         gene_sd_floor = 0.10,
                                         score_sd_floor = 0.10,
                                         min_matched_genes = 3L) {
  prep <- prepare_dataset_expression(dataset_id)
  symbol_map <- load_gene_symbol_map() %>%
    filter(gene_id %in% rownames(prep$x)) %>%
    mutate(gene_symbol_lower = tolower(gene_symbol))
  ctrl_ids <- prep$design$sample_id[prep$design$condition_simple == "CONTROL"]

  sig_qc <- list()
  sig_scores <- list()
  for (sig in unique(signature_defs$signature_id)) {
    sig_def <- signature_defs %>% filter(signature_id == sig)
    matched <- symbol_map %>%
      inner_join(sig_def %>% select(signature_id, gene_symbol_lower),
                 by = "gene_symbol_lower") %>%
      distinct(gene_id, .keep_all = TRUE)
    requested_symbols <- sig_def$gene_symbol
    matched_symbols <- unique(matched$gene_symbol)
    missing_symbols <- setdiff(requested_symbols, matched_symbols)
    n_requested <- length(unique(requested_symbols))
    n_matched <- nrow(matched)
    signature_label <- sig_def$signature_label[[1]]

    sig_qc[[sig]] <- tibble(
      dataset_id = dataset_id,
      signature_id = sig,
      signature_label = signature_label,
      n_requested_genes = n_requested,
      n_matched_genes = n_matched,
      coverage_fraction = n_matched / n_requested,
      matched_gene_symbols = paste(sort(matched_symbols), collapse = ";"),
      missing_gene_symbols = paste(sort(missing_symbols), collapse = ";"),
      counts_path = prep$counts_path,
      design_path = prep$design_path,
      status = ifelse(n_matched >= min_matched_genes, "ok",
                      "failed_too_few_matched_genes")
    )

    if (n_matched < min_matched_genes) next
    xsub <- prep$x[matched$gene_id, , drop = FALSE]
    ctrl_mat <- xsub[, ctrl_ids, drop = FALSE]
    mu <- rowMeans(ctrl_mat, na.rm = TRUE)
    sdv <- apply(ctrl_mat, 1, sd, na.rm = TRUE)
    sdv_floor <- pmax(sdv, gene_sd_floor)
    z <- sweep(xsub, 1, mu, "-")
    z <- sweep(z, 1, sdv_floor, "/")
    raw <- colMeans(z, na.rm = TRUE)
    ctrl_raw <- raw[ctrl_ids]
    raw_sd <- stats::sd(ctrl_raw, na.rm = TRUE)
    raw_sd_floor <- max(raw_sd, score_sd_floor, na.rm = TRUE)
    if (!is.finite(raw_sd_floor) || raw_sd_floor <= 0) raw_sd_floor <- score_sd_floor
    score_z <- (raw - mean(ctrl_raw, na.rm = TRUE)) / raw_sd_floor

    sig_scores[[sig]] <- prep$design %>%
      transmute(
        dataset_id = dataset_id,
        sample_id = sample_id,
        condition_simple = condition_simple,
        signature_id = sig,
        signature_label = signature_label,
        signature_raw = as.numeric(raw[sample_id]),
        signature_z = as.numeric(score_z[sample_id]),
        n_requested_genes = n_requested,
        n_matched_genes = n_matched,
        coverage_fraction = n_matched / n_requested
      )
  }

  list(scores = bind_rows(sig_scores), qc = bind_rows(sig_qc))
}

evaluate_signature_for_split <- function(score_tbl, row, signature_id) {
  split <- read_split_table(row$split_path)
  one_score <- score_tbl %>% filter(signature_id == !!signature_id)
  merged <- split %>%
    select(any_of(c("sample_id", "condition_simple", "tissue", "time_h",
                    "control_label", "contrast_label"))) %>%
    inner_join(one_score, by = "sample_id", suffix = c("_split", "")) %>%
    mutate(
      condition_simple = ifelse(!is.na(condition_simple_split),
                                condition_simple_split, condition_simple),
      label = ifelse(condition_simple == "DELIVERY", 1L, 0L)
    )

  ctrl <- merged$signature_z[merged$condition_simple == "CONTROL"]
  del <- merged$signature_z[merged$condition_simple == "DELIVERY"]
  pass <- length(ctrl) > 0 && length(del) > 0
  tibble(
    gse_id = row$gse_id,
    dataset_id = row$gse_id,
    split_id = row$split_id,
    split_path = row$split_path,
    manuscript_role = row$manuscript_role,
    manuscript_interpretation_group = row$manuscript_interpretation_group,
    manuscript_interpretation_label = row$manuscript_interpretation_label,
    contrast_label = row$contrast_label,
    control_label = row$control_label,
    tissue = row$tissue,
    time_h = row$time_h,
    signature_id = signature_id,
    signature_label = unique(one_score$signature_label)[[1]],
    pass = pass,
    fail_reason = ifelse(pass, NA_character_, "missing_control_or_delivery"),
    n_scored_samples = nrow(merged),
    n_controls = sum(merged$condition_simple == "CONTROL"),
    n_delivery = sum(merged$condition_simple == "DELIVERY"),
    n_requested_genes = unique(one_score$n_requested_genes)[[1]],
    n_matched_genes = unique(one_score$n_matched_genes)[[1]],
    coverage_fraction = unique(one_score$coverage_fraction)[[1]],
    control_mean_signature_z = safe_mean(ctrl),
    delivery_mean_signature_z = safe_mean(del),
    delta_mean_signature_z = safe_mean(del) - safe_mean(ctrl),
    control_median_signature_z = safe_median(ctrl),
    delivery_median_signature_z = safe_median(del),
    delta_median_signature_z = safe_median(del) - safe_median(ctrl),
    cohens_d_signature_z = cohens_d_two_group(del, ctrl),
    t_test_p_signature_z = ifelse(
      length(ctrl) >= 2 && length(del) >= 2,
      tryCatch(stats::t.test(del, ctrl)$p.value, error = function(e) NA_real_),
      NA_real_
    ),
    auc_signature_z = compute_auc_safe(merged$signature_z, merged$label)
  )
}

main <- function() {
  log_extra("Starting baseline signature benchmarking.")
  eval_tbl <- load_eval_with_roles(required = TRUE) %>%
    filter(pass %in% TRUE)
  if (nrow(eval_tbl) == 0) stop("No passing Step 09 contrasts were available.")
  signature_defs <- make_signature_definitions()

  write_tsv_safe(signature_defs %>%
                   select(signature_id, signature_label, gene_symbol),
                 file.path(EXTRA_RESULTS_DIR, "baseline_signature_requested_genes.tsv"))

  datasets <- sort(unique(eval_tbl$gse_id))
  scored <- map(datasets, function(ds) {
    log_extra("Scoring baseline signatures for ", ds, ".")
    tryCatch(
      score_signatures_for_dataset(ds, signature_defs),
      error = function(e) {
        log_extra("WARN: baseline scoring failed for ", ds, ": ",
                  conditionMessage(e))
        list(scores = tibble(), qc = tibble(
          dataset_id = ds,
          signature_id = NA_character_,
          signature_label = NA_character_,
          n_requested_genes = NA_integer_,
          n_matched_genes = NA_integer_,
          coverage_fraction = NA_real_,
          matched_gene_symbols = NA_character_,
          missing_gene_symbols = NA_character_,
          counts_path = find_counts_file(ds),
          design_path = find_scoring_design(ds),
          status = paste0("failed: ", conditionMessage(e))
        ))
      }
    )
  })
  names(scored) <- datasets
  sample_scores <- bind_rows(map(scored, "scores"))
  qc_tbl <- bind_rows(map(scored, "qc"))

  matched_defs <- qc_tbl %>%
    filter(status == "ok") %>%
    select(dataset_id, signature_id, signature_label, n_requested_genes,
           n_matched_genes, coverage_fraction, matched_gene_symbols,
           missing_gene_symbols) %>%
    arrange(dataset_id, signature_id)

  eval_rows <- pmap_dfr(
    list(eval_tbl$gse_id, seq_len(nrow(eval_tbl))),
    function(ds, idx) {
      row <- eval_tbl[idx, ]
      scores_ds <- sample_scores %>% filter(dataset_id == ds)
      if (nrow(scores_ds) == 0) return(tibble())
      map_dfr(unique(scores_ds$signature_id), function(sig) {
        tryCatch(
          evaluate_signature_for_split(scores_ds, row, sig),
          error = function(e) {
            tibble(
              gse_id = row$gse_id,
              dataset_id = row$gse_id,
              split_id = row$split_id,
              split_path = row$split_path,
              manuscript_role = row$manuscript_role,
              manuscript_interpretation_group =
                row$manuscript_interpretation_group,
              manuscript_interpretation_label =
                row$manuscript_interpretation_label,
              signature_id = sig,
              pass = FALSE,
              fail_reason = conditionMessage(e)
            )
          }
        )
      })
    }
  )

  imrs_contrast <- eval_tbl %>%
    transmute(
      gse_id, split_id, split_path, manuscript_role,
      manuscript_interpretation_group, manuscript_interpretation_label,
      contrast_label, control_label, tissue, time_h,
      score_id = "IMRS",
      score_label = "IMRS",
      delta_score = delta_mean_imrs_z,
      cohens_d = cohens_d_imrs_z,
      auc_secondary = auc_imrs_z,
      n_controls, n_delivery,
      pass = TRUE
    )
  baseline_contrast <- eval_rows %>%
    filter(pass %in% TRUE) %>%
    transmute(
      gse_id, split_id, split_path, manuscript_role,
      manuscript_interpretation_group, manuscript_interpretation_label,
      contrast_label, control_label, tissue, time_h,
      score_id = signature_id,
      score_label = signature_label,
      delta_score = delta_mean_signature_z,
      cohens_d = cohens_d_signature_z,
      auc_secondary = auc_signature_z,
      n_controls, n_delivery,
      pass = TRUE
    )
  contrast_long <- bind_rows(imrs_contrast, baseline_contrast) %>%
    mutate(
      manuscript_interpretation_group = factor(
        as.character(manuscript_interpretation_group),
        levels = INTERPRETATION_ORDER
      ),
      score_label = factor(
        score_label,
        levels = c("IMRS", unique(signature_defs$signature_label))
      )
    )

  summary_by_group <- contrast_long %>%
    filter(pass %in% TRUE) %>%
    group_by(score_id, score_label, manuscript_interpretation_group,
             manuscript_interpretation_label) %>%
    summarise(
      n_contrasts = n(),
      mean_delta_score = safe_mean(delta_score),
      median_delta_score = safe_median(delta_score),
      proportion_positive_delta = safe_prop(delta_score > 0),
      mean_auc_secondary = safe_mean(auc_secondary),
      median_auc_secondary = safe_median(auc_secondary),
      .groups = "drop"
    ) %>%
    arrange(score_label,
            match(as.character(manuscript_interpretation_group),
                  INTERPRETATION_ORDER))

  paired_contrast <- baseline_contrast %>%
    select(gse_id, split_id, signature_id = score_id,
           signature_label = score_label,
           baseline_delta_score = delta_score,
           baseline_auc_secondary = auc_secondary) %>%
    left_join(imrs_contrast %>%
                select(gse_id, split_id, imrs_delta_score = delta_score,
                       imrs_auc_secondary = auc_secondary),
              by = c("gse_id", "split_id"))

  imrs_sample <- load_sample_with_roles(required = TRUE) %>%
    transmute(gse_id, split_id, sample_id, condition_simple, imrs_z,
              manuscript_interpretation_group, manuscript_interpretation_label)
  baseline_sample_split <- imrs_sample %>%
    select(gse_id, split_id, sample_id, condition_simple,
           manuscript_interpretation_group, manuscript_interpretation_label,
           imrs_z) %>%
    inner_join(sample_scores %>%
                 select(dataset_id, sample_id, signature_id, signature_label,
                        signature_z),
               by = c("gse_id" = "dataset_id", "sample_id"),
               relationship = "many-to-many")
  correlation_summary <- baseline_sample_split %>%
    group_by(signature_id, signature_label) %>%
    summarise(
      n_split_sample_records = sum(is.finite(imrs_z) & is.finite(signature_z)),
      pearson_r = suppressWarnings(stats::cor(imrs_z, signature_z,
                                              use = "complete.obs",
                                              method = "pearson")),
      spearman_r = suppressWarnings(stats::cor(imrs_z, signature_z,
                                               use = "complete.obs",
                                               method = "spearman")),
      .groups = "drop"
    )

  write_tsv_safe(qc_tbl, file.path(EXTRA_RESULTS_DIR,
                                   "baseline_signature_dataset_qc.tsv"))
  write_tsv_safe(matched_defs, file.path(EXTRA_RESULTS_DIR,
                                         "baseline_signature_gene_sets.tsv"))
  write_tsv_safe(sample_scores, file.path(EXTRA_RESULTS_DIR,
                                          "baseline_signature_scores_sample_level.tsv"))
  write_tsv_safe(eval_rows, file.path(EXTRA_RESULTS_DIR,
                                      "baseline_signature_contrast_eval.tsv"))
  write_csv_safe(eval_rows, file.path(EXTRA_RESULTS_DIR,
                                      "baseline_signature_contrast_eval.csv"))
  write_tsv_safe(contrast_long, file.path(EXTRA_RESULTS_DIR,
                                          "baseline_signature_contrast_long.tsv"))
  write_tsv_safe(summary_by_group, file.path(EXTRA_RESULTS_DIR,
                                             "baseline_signature_summary_by_group.tsv"))
  write_tsv_safe(paired_contrast, file.path(EXTRA_RESULTS_DIR,
                                            "baseline_signature_paired_contrast_comparison.tsv"))
  write_tsv_safe(correlation_summary, file.path(EXTRA_RESULTS_DIR,
                                                "baseline_signature_correlation_with_imrs.tsv"))

  p_delta <- ggplot(contrast_long %>% filter(pass %in% TRUE),
                    aes(x = manuscript_interpretation_group,
                        y = delta_score,
                        fill = manuscript_interpretation_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72) +
    geom_jitter(width = 0.12, height = 0, size = 1.35, alpha = 0.62) +
    facet_wrap(~ score_label, scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = INTERPRETATION_LABELS, drop = FALSE) +
    scale_fill_manual(values = INTERPRETATION_PALETTE, drop = FALSE,
                      labels = INTERPRETATION_LABELS) +
    labs(
      title = "Delta Score by Manuscript Interpretation Group",
      subtitle = "Baseline signatures use unweighted control-standardized gene z-score averages.",
      x = NULL,
      y = "Delta score (DELIVERY - CONTROL)",
      fill = "Manuscript group"
    ) +
    theme_extra() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_plot_pair(p_delta, "baseline_delta_by_manuscript_group",
                 width = 9, height = 6.4)

  p_paired <- ggplot(paired_contrast,
                     aes(x = imrs_delta_score, y = baseline_delta_score,
                         color = signature_label)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_point(alpha = 0.72, size = 1.7) +
    facet_wrap(~ signature_label, scales = "free_y") +
    labs(
      title = "Paired Contrast-Level Delta Comparison",
      subtitle = "Each point is the same split contrast scored by IMRS and by a simple baseline signature.",
      x = "IMRS Delta",
      y = "Baseline signature Delta",
      color = "Baseline signature"
    ) +
    theme_extra()
  save_plot_pair(p_paired, "baseline_paired_contrast_comparison",
                 width = 8.5, height = 5.5)

  direction_tbl <- summary_by_group %>%
    filter(!is.na(proportion_positive_delta))
  p_direction <- ggplot(direction_tbl,
                        aes(x = manuscript_interpretation_group,
                            y = proportion_positive_delta,
                            fill = score_label)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_x_discrete(labels = INTERPRETATION_LABELS, drop = FALSE) +
    labs(
      title = "Directionality Summary",
      subtitle = "Fraction of split contrasts with positive DELIVERY - CONTROL delta.",
      x = NULL,
      y = "Positive-delta contrasts",
      fill = "Score"
    ) +
    theme_extra() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_plot_pair(p_direction, "baseline_directionality_summary",
                 width = 8.8, height = 5.2)

  p_corr <- ggplot(baseline_sample_split,
                   aes(x = signature_z, y = imrs_z,
                       color = manuscript_interpretation_group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.45) +
    facet_wrap(~ signature_label, scales = "free") +
    scale_color_manual(values = INTERPRETATION_PALETTE, drop = FALSE,
                       labels = INTERPRETATION_LABELS) +
    labs(
      title = "Sample-Level Correlation of IMRS with Baseline Signatures",
      subtitle = "Records are split-level samples, so shared controls may appear in more than one contrast.",
      x = "Baseline signature z-score",
      y = "IMRSz",
      color = "Manuscript group"
    ) +
    theme_extra()
  save_plot_pair(p_corr, "baseline_correlation_with_imrs",
                 width = 8.6, height = 5.6)

  log_extra("Completed baseline signature benchmarking. Contrast rows=",
            nrow(eval_rows), ".")
  invisible(summary_by_group)
}

main()
