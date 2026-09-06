#!/usr/bin/env Rscript

# Label permutation/null analysis for frozen IMRS scores.
# Important: labels are permuted within each existing split contrast while
# imrs_z values are kept fixed. No weights are retrained or re-estimated.

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

this_script <- imrs_detect_script_path_local("01_label_permutation_null.R")
if (is.na(this_script)) stop("Could not identify this script path. Use RStudio Source/Run or Rscript.", call. = FALSE)
source(file.path(dirname(this_script), "00_publication_extra_config.R"))

set.seed(20260427)

args <- commandArgs(trailingOnly = TRUE)
n_perm <- if (length(args) >= 1) suppressWarnings(as.integer(args[[1]])) else 1000L
if (!is.finite(n_perm) || n_perm < 100L) n_perm <- 1000L

permutation_delta <- function(scores, labels, n_delivery, n_perm) {
  n <- length(scores)
  replicate(n_perm, {
    delivery_idx <- sample.int(n, size = n_delivery, replace = FALSE)
    mean(scores[delivery_idx], na.rm = TRUE) -
      mean(scores[-delivery_idx], na.rm = TRUE)
  })
}

summarise_one_split <- function(df, n_perm) {
  df <- df %>%
    filter(condition_simple %in% c("CONTROL", "DELIVERY"),
           is.finite(imrs_z)) %>%
    distinct(sample_id, .keep_all = TRUE)
  n_control <- sum(df$condition_simple == "CONTROL")
  n_delivery <- sum(df$condition_simple == "DELIVERY")
  if (n_control < 1 || n_delivery < 1 || nrow(df) < 3) {
    return(list(summary = tibble(
      permutation_status = "failed",
      fail_reason = "insufficient_control_or_delivery"
    ), null = tibble()))
  }

  scores <- df$imrs_z
  labels <- df$condition_simple
  observed <- mean(scores[labels == "DELIVERY"], na.rm = TRUE) -
    mean(scores[labels == "CONTROL"], na.rm = TRUE)
  null_delta <- permutation_delta(scores, labels, n_delivery, n_perm)
  qs <- safe_quantile(null_delta, c(0.025, 0.05, 0.5, 0.95, 0.975))
  null_sd <- safe_sd(null_delta)
  empirical_p_two_sided <- (1 + sum(abs(null_delta) >= abs(observed), na.rm = TRUE)) /
    (n_perm + 1)
  empirical_p_greater <- (1 + sum(null_delta >= observed, na.rm = TRUE)) /
    (n_perm + 1)
  empirical_p_less <- (1 + sum(null_delta <= observed, na.rm = TRUE)) /
    (n_perm + 1)

  list(
    summary = tibble(
      permutation_status = "ok",
      fail_reason = NA_character_,
      n_permutations = n_perm,
      n_control = n_control,
      n_delivery = n_delivery,
      observed_delta_mean_imrs_z = observed,
      null_mean_delta = safe_mean(null_delta),
      null_sd_delta = null_sd,
      null_q025_delta = qs[[1]],
      null_q05_delta = qs[[2]],
      null_q50_delta = qs[[3]],
      null_q95_delta = qs[[4]],
      null_q975_delta = qs[[5]],
      null_z = ifelse(is.finite(null_sd) && null_sd > 0,
                      (observed - safe_mean(null_delta)) / null_sd,
                      NA_real_),
      empirical_p_two_sided = empirical_p_two_sided,
      empirical_p_greater = empirical_p_greater,
      empirical_p_less = empirical_p_less
    ),
    null = tibble(
      permutation_index = seq_len(n_perm),
      permuted_delta_mean_imrs_z = null_delta
    )
  )
}

main <- function() {
  log_extra("Starting label permutation/null analysis with ", n_perm,
            " permutations per valid contrast.")

  eval_tbl <- load_eval_with_roles(required = TRUE) %>%
    filter(pass %in% TRUE)
  sample_tbl <- load_sample_with_roles(required = TRUE) %>%
    semi_join(eval_tbl %>% select(gse_id, split_id), by = c("gse_id", "split_id"))
  if (nrow(sample_tbl) == 0) stop("No valid split sample rows were available.")

  split_keys <- sample_tbl %>%
    distinct(gse_id, split_id, split_path, manuscript_role,
             manuscript_interpretation_group, manuscript_interpretation_label,
             contrast_label, control_label, tissue, time_h) %>%
    arrange(gse_id, split_id)

  results <- vector("list", nrow(split_keys))
  null_rows <- vector("list", nrow(split_keys))
  for (i in seq_len(nrow(split_keys))) {
    key <- split_keys[i, ]
    df <- sample_tbl %>%
      filter(gse_id == key$gse_id, split_id == key$split_id)
    out <- summarise_one_split(df, n_perm)
    results[[i]] <- bind_cols(key, out$summary)
    if (nrow(out$null) > 0) {
      null_rows[[i]] <- bind_cols(
        key %>% select(gse_id, split_id, manuscript_interpretation_group,
                       manuscript_interpretation_label),
        out$null
      )
    } else {
      null_rows[[i]] <- tibble()
    }
  }

  summary_tbl <- bind_rows(results) %>%
    mutate(
      empirical_p_two_sided_fdr = p.adjust(empirical_p_two_sided,
                                           method = "BH"),
      empirical_p_greater_fdr = p.adjust(empirical_p_greater, method = "BH"),
      empirical_p_less_fdr = p.adjust(empirical_p_less, method = "BH"),
      observed_outside_95pct_null =
        observed_delta_mean_imrs_z < null_q025_delta |
        observed_delta_mean_imrs_z > null_q975_delta
    ) %>%
    arrange(manuscript_interpretation_group, gse_id, time_h, split_id)
  null_tbl <- bind_rows(null_rows)

  group_summary <- summary_tbl %>%
    filter(permutation_status == "ok") %>%
    group_by(manuscript_interpretation_group, manuscript_interpretation_label) %>%
    summarise(
      n_contrasts = n(),
      mean_observed_delta = safe_mean(observed_delta_mean_imrs_z),
      median_observed_delta = safe_median(observed_delta_mean_imrs_z),
      proportion_positive_delta = safe_prop(observed_delta_mean_imrs_z > 0),
      median_empirical_p_greater = safe_median(empirical_p_greater),
      median_empirical_p_two_sided = safe_median(empirical_p_two_sided),
      proportion_fdr_greater_below_0_05 =
        safe_prop(empirical_p_greater_fdr < 0.05),
      proportion_outside_95pct_null = safe_prop(observed_outside_95pct_null),
      .groups = "drop"
    ) %>%
    arrange(match(as.character(manuscript_interpretation_group),
                  INTERPRETATION_ORDER))

  write_tsv_safe(summary_tbl, file.path(EXTRA_RESULTS_DIR,
                                        "label_permutation_null_summary.tsv"))
  write_csv_safe(summary_tbl, file.path(EXTRA_RESULTS_DIR,
                                        "label_permutation_null_summary.csv"))
  write_tsv_safe(null_tbl, file.path(EXTRA_RESULTS_DIR,
                                     "label_permutation_null_long.tsv"))
  write_tsv_safe(group_summary, file.path(EXTRA_RESULTS_DIR,
                                          "label_permutation_group_summary.tsv"))

  plot_tbl <- summary_tbl %>%
    filter(permutation_status == "ok") %>%
    mutate(
      split_short = paste(gse_id, row_number(), sep = "_"),
      order_index = row_number()
    ) %>%
    arrange(observed_delta_mean_imrs_z) %>%
    mutate(order_index = row_number())

  p_observed <- ggplot(plot_tbl, aes(x = order_index)) +
    geom_linerange(aes(ymin = null_q025_delta, ymax = null_q975_delta),
                   color = "#BDBDBD", linewidth = 0.7) +
    geom_linerange(aes(ymin = null_q05_delta, ymax = null_q95_delta),
                   color = "#777777", linewidth = 0.35) +
    geom_point(aes(y = observed_delta_mean_imrs_z,
                   color = manuscript_interpretation_group), size = 1.7,
               alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    scale_color_manual(values = INTERPRETATION_PALETTE, drop = FALSE,
                       labels = INTERPRETATION_LABELS) +
    labs(
      title = "Observed Delta IMRSz Versus Label-Permutation Null",
      subtitle = "Gray ranges show per-contrast null intervals; points show frozen-score observed deltas.",
      x = "Split contrasts ordered by observed Delta IMRSz",
      y = "Delta IMRSz",
      color = "Manuscript group"
    ) +
    theme_extra()
  save_plot_pair(p_observed, "label_permutation_observed_vs_null",
                 width = 8.4, height = 5.1)

  p_pvals <- ggplot(summary_tbl %>% filter(permutation_status == "ok"),
                    aes(x = empirical_p_greater)) +
    geom_histogram(binwidth = 0.05, boundary = 0, color = "white",
                   fill = "#4E79A7") +
    facet_wrap(~ manuscript_interpretation_label, ncol = 2) +
    labs(
      title = "Empirical Activation P-Value Distribution",
      subtitle = "One-sided label-permutation p-values for observed positive delivery-associated Delta IMRSz.",
      x = "Empirical p-value, greater-than null",
      y = "Number of split contrasts"
    ) +
    theme_extra()
  save_plot_pair(p_pvals, "label_permutation_empirical_pvalue_distribution",
                 width = 8, height = 5.4)

  p_groups <- ggplot(summary_tbl %>% filter(permutation_status == "ok"),
                     aes(x = manuscript_interpretation_group,
                         y = observed_delta_mean_imrs_z,
                         fill = manuscript_interpretation_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72) +
    geom_jitter(width = 0.13, height = 0, size = 1.7, alpha = 0.75) +
    scale_x_discrete(labels = INTERPRETATION_LABELS, drop = FALSE) +
    scale_fill_manual(values = INTERPRETATION_PALETTE, drop = FALSE,
                      labels = INTERPRETATION_LABELS) +
    labs(
      title = "Permutation Results by Manuscript Interpretation Group",
      subtitle = "Primary acute validation is kept separate from extended and secondary support.",
      x = NULL,
      y = "Observed Delta IMRSz",
      fill = "Manuscript group"
    ) +
    theme_extra() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_plot_pair(p_groups, "label_permutation_primary_vs_extended_secondary",
                 width = 8.2, height = 5.2)

  log_extra("Completed label permutation/null analysis. Contrasts=", nrow(summary_tbl),
            "; null rows=", nrow(null_tbl), ".")
  invisible(summary_tbl)
}

main()
