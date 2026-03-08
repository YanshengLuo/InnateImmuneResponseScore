#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 09 — Calibration Evaluation by SPLIT CONTRAST
#
# Corrected version:
#   1) output dataset_id is based on DELIVERY sample_id from split file
#      (usually SRR...)
#   2) time_h is taken from:
#        a) H= in split filename, or
#        b) DELIVERY-only time_h
#      and never averaged with control 0h samples
#   3) split files are searched recursively
#   4) only true contrast files containing "__VS=" are used
#   5) plotting guards added to avoid crashes on sparse datasets
#
# Inputs:
#   1) Step 08 score files:
#      <project_root>/05_score/transfer/scores/*__imrs_scores.tsv
#
#   2) Verified split design files:
#      <project_root>/00_metadata/verified_metadata/splited/<DATASET or DATASET_design>/**/*.tsv
#
# Outputs:
#   Tables:
#     <project_root>/05_score/transfer/eval/step09_split_eval.tsv
#     <project_root>/05_score/transfer/eval/step09_split_summary.tsv
#     <project_root>/05_score/transfer/eval/step09_split_sample_level.tsv
#
#   Figures:
#     <project_root>/06_figures/step09_calibration/
#
# Notes:
#   - Evaluation only
#   - No retraining
#   - No weight updates
#   - Primary metric: delta_mean_imrs_z = mean(DELIVERY) - mean(CONTROL)
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

# Edit this to the datasets you want to evaluate
calibration_dataset_ids <- c(
"GSE262515"
)

do_auc <- TRUE
write_roc <- TRUE

scores_dir <- file.path(project_root, "05_score", "transfer", "scores")
split_root <- file.path(project_root, "00_metadata", "verified_metadata", "splited")
eval_dir   <- file.path(project_root, "05_score", "transfer", "eval")
roc_dir    <- file.path(eval_dir, "roc")
fig_dir    <- file.path(project_root, "06_figures", "step09_calibration")

dir.create(eval_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(roc_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(scores_dir)) stop("Missing scores_dir: ", scores_dir)
if (!dir.exists(split_root)) stop("Missing split_root: ", split_root)

# -------------------------
# HELPERS
# -------------------------
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

collapse_unique_chr <- function(x, sep = ";") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

cohens_d <- function(x_treat, x_ctrl) {
  x_treat <- x_treat[is.finite(x_treat)]
  x_ctrl  <- x_ctrl[is.finite(x_ctrl)]

  if (length(x_treat) < 2 || length(x_ctrl) < 2) return(NA_real_)

  m1 <- mean(x_treat)
  m0 <- mean(x_ctrl)
  s1 <- sd(x_treat)
  s0 <- sd(x_ctrl)

  sp <- sqrt(
    (((length(x_treat) - 1) * s1^2) + ((length(x_ctrl) - 1) * s0^2)) /
      (length(x_treat) + length(x_ctrl) - 2)
  )

  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m0) / sp
}

extract_gse_id_from_scorefile <- function(path) {
  b <- basename(path)
  sub("__.*$", "", b)
}

extract_score_id <- function(path) {
  sub("__imrs_scores\\.tsv$", "", basename(path))
}

read_split_table <- function(path) {
  out <- tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) stop("Could not read split file: ", path)
  out
}

extract_H_from_split_filename <- function(path) {
  b <- sub("\\.tsv$", "", basename(path))
  m <- str_match(b, "__H=([0-9.]+)(?:__|$)")
  val <- suppressWarnings(as.numeric(m[, 2]))
  if (length(val) == 0 || is.na(val)) return(NA_real_)
  val
}

compute_eval_time_h <- function(split_path, merged_df) {
  h_from_name <- extract_H_from_split_filename(split_path)
  if (is.finite(h_from_name)) return(h_from_name)

  if ("time_h" %in% names(merged_df)) {
    del_times <- merged_df$time_h[merged_df$condition_simple == "DELIVERY"]
    del_times <- suppressWarnings(as.numeric(del_times))
    del_times <- del_times[is.finite(del_times)]
    if (length(del_times) > 0) return(mean(unique(del_times)))
  }

  NA_real_
}

compute_auc <- function(df) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    return(list(auc = NA_real_, roc_tbl = NULL, auc_note = "pROC_not_installed"))
  }

  tmp <- df %>%
    filter(is.finite(imrs_z), !is.na(label))

  if (nrow(tmp) == 0 || length(unique(tmp$label)) < 2) {
    return(list(auc = NA_real_, roc_tbl = NULL, auc_note = "one_class_only"))
  }

  roc_obj <- pROC::roc(
    response = tmp$label,
    predictor = tmp$imrs_z,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_tbl <- tibble(
    threshold   = roc_obj$thresholds,
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )

  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    roc_tbl = roc_tbl,
    auc_note = "ok"
  )
}

choose_best_score_file <- function(x) {
  if (length(x) == 1) return(x)
  preferred <- x[grepl("__counts__imrs_scores\\.tsv$", basename(x))]
  if (length(preferred) >= 1) return(preferred[1])
  x[1]
}

save_plot_both <- function(plot_obj, out_stub, width = 8, height = 6) {
  ggsave(
    filename = paste0(out_stub, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    filename = paste0(out_stub, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# -------------------------
# INDEX STEP 08 SCORE FILES
# -------------------------
all_score_files <- list.files(
  scores_dir,
  pattern = "__imrs_scores\\.tsv$",
  full.names = TRUE
)

if (length(all_score_files) == 0) {
  stop("No Step 08 score files found under: ", scores_dir)
}

score_index <- tibble(
  scores_path = all_score_files,
  gse_id = vapply(all_score_files, extract_gse_id_from_scorefile, character(1)),
  score_id = vapply(all_score_files, extract_score_id, character(1))
) %>%
  filter(gse_id %in% calibration_dataset_ids)

if (nrow(score_index) == 0) {
  stop("No Step 08 scores found for calibration datasets: ",
       paste(calibration_dataset_ids, collapse = ", "))
}

score_map <- score_index %>%
  group_by(gse_id) %>%
  summarise(
    scores_path = choose_best_score_file(scores_path),
    score_id = choose_best_score_file(score_id),
    .groups = "drop"
  )

# -------------------------
# MAIN EVALUATION
# -------------------------
eval_rows <- list()
sample_level_rows <- list()

for (ds in calibration_dataset_ids) {

  split_dir_candidates <- c(
    file.path(split_root, ds),
    file.path(split_root, paste0(ds, "_design"))
  )

  existing_split_dirs <- split_dir_candidates[dir.exists(split_dir_candidates)]

  if (length(existing_split_dirs) == 0) {
    message("Skip dataset (missing split dir): ", ds)
    message("  tried: ", paste(split_dir_candidates, collapse = " | "))
    next
  }

  ds_split_dir <- existing_split_dirs[1]

  score_row <- score_map %>% filter(gse_id == ds)
  if (nrow(score_row) == 0) {
    message("Skip dataset (no Step 08 score file): ", ds)
    next
  }

  score_path <- score_row$scores_path[[1]]

  scores_df <- read_tsv(score_path, show_col_types = FALSE, progress = FALSE)

  required_score_cols <- c("sample_id", "condition_simple", "imrs_raw", "imrs_z", "coverage")
  miss_score <- setdiff(required_score_cols, names(scores_df))
  if (length(miss_score) > 0) {
    stop("Step 08 score file missing columns for ", ds, ": ",
         paste(miss_score, collapse = ", "))
  }

  scores_df <- scores_df %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple))),
      imrs_raw = suppressWarnings(as.numeric(imrs_raw)),
      imrs_z   = suppressWarnings(as.numeric(imrs_z)),
      coverage = suppressWarnings(as.numeric(coverage))
    )

  split_files <- list.files(
    ds_split_dir,
    pattern = "\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )

  # keep only true contrast split files
  split_files <- split_files[grepl("__VS=", basename(split_files))]

  if (length(split_files) == 0) {
    message("Skip dataset (no split contrast TSVs): ", ds)
    next
  }

  message("Dataset: ", ds, " | split files: ", length(split_files))
  # message("First few split files:")
  # print(head(split_files))

  for (sf in split_files) {
    split_id <- sub("\\.tsv$", "", basename(sf))

    split_df <- tryCatch(read_split_table(sf), error = function(e) NULL)
    if (is.null(split_df)) {
      eval_rows[[length(eval_rows) + 1]] <- tibble(
        gse_id = ds,
        dataset_id = NA_character_,
        split_id = split_id,
        split_path = sf,
        pass = FALSE,
        fail_reason = "Could not read split file"
      )
      next
    }

    required_split_cols <- c(
      "sample_id", "condition_simple", "contrast_label",
      "control_label", "tissue"
    )
    miss_split <- setdiff(required_split_cols, names(split_df))
    if (length(miss_split) > 0) {
      eval_rows[[length(eval_rows) + 1]] <- tibble(
        gse_id = ds,
        dataset_id = NA_character_,
        split_id = split_id,
        split_path = sf,
        pass = FALSE,
        fail_reason = paste0("Split file missing columns: ", paste(miss_split, collapse = ", "))
      )
      next
    }

    if (!("time_h" %in% names(split_df))) split_df$time_h <- NA_real_

    split_df <- split_df %>%
      mutate(
        sample_id = trimws(as.character(sample_id)),
        condition_simple = toupper(trimws(as.character(condition_simple))),
        contrast_label = as.character(contrast_label),
        control_label = as.character(control_label),
        tissue = as.character(tissue),
        time_h = suppressWarnings(as.numeric(time_h))
      )

    merged <- split_df %>%
      inner_join(scores_df, by = c("sample_id", "condition_simple"), suffix = c("_split", "_score"))

    n_split_samples <- nrow(split_df)
    n_scored_samples <- nrow(merged)

    ctrl <- merged %>% filter(condition_simple == "CONTROL")
    del  <- merged %>% filter(condition_simple == "DELIVERY")

    n_ctrl <- nrow(ctrl)
    n_del  <- nrow(del)

    dataset_id_out <- collapse_unique_chr(del$sample_id)

    if (n_scored_samples == 0) {
      eval_rows[[length(eval_rows) + 1]] <- tibble(
        gse_id = ds,
        dataset_id = dataset_id_out,
        split_id = split_id,
        split_path = sf,
        pass = FALSE,
        fail_reason = "No overlap between split samples and scored samples",
        n_split_samples = n_split_samples,
        n_scored_samples = 0
      )
      next
    }

    if (n_ctrl < 1 || n_del < 1) {
      eval_rows[[length(eval_rows) + 1]] <- tibble(
        gse_id = ds,
        dataset_id = dataset_id_out,
        split_id = split_id,
        split_path = sf,
        pass = FALSE,
        fail_reason = "Need both CONTROL and DELIVERY in split",
        n_split_samples = n_split_samples,
        n_scored_samples = n_scored_samples,
        n_controls = n_ctrl,
        n_delivery = n_del
      )
      next
    }

    ctrl_mean   <- safe_mean(ctrl$imrs_z)
    del_mean    <- safe_mean(del$imrs_z)
    ctrl_median <- safe_median(ctrl$imrs_z)
    del_median  <- safe_median(del$imrs_z)

    delta_mean   <- del_mean - ctrl_mean
    delta_median <- del_median - ctrl_median
    d_val        <- cohens_d(del$imrs_z, ctrl$imrs_z)

    t_p <- tryCatch(
      t.test(del$imrs_z, ctrl$imrs_z)$p.value,
      error = function(e) NA_real_
    )

    auc_val <- NA_real_
    auc_note <- NA_character_

    if (isTRUE(do_auc)) {
      auc_input <- merged %>%
        mutate(label = ifelse(condition_simple == "DELIVERY", 1L, 0L))

      auc_out <- compute_auc(auc_input)
      auc_val <- auc_out$auc
      auc_note <- auc_out$auc_note

      if (isTRUE(write_roc) && !is.null(auc_out$roc_tbl)) {
        write_tsv(
          auc_out$roc_tbl,
          file.path(roc_dir, paste0(split_id, "__roc.tsv"))
        )
      }
    }

    eval_time_h <- compute_eval_time_h(sf, merged)

    eval_rows[[length(eval_rows) + 1]] <- tibble(
      gse_id = ds,
      dataset_id = dataset_id_out,
      split_id = split_id,
      split_path = sf,
      pass = TRUE,
      fail_reason = NA_character_,

      contrast_label = collapse_unique_chr(merged$contrast_label),
      control_label = collapse_unique_chr(merged$control_label),
      tissue = collapse_unique_chr(merged$tissue),
      time_h = eval_time_h,

      delivery_sample_ids = collapse_unique_chr(del$sample_id),
      control_sample_ids  = collapse_unique_chr(ctrl$sample_id),

      n_split_samples = n_split_samples,
      n_scored_samples = n_scored_samples,
      n_controls = n_ctrl,
      n_delivery = n_del,
      mean_coverage = safe_mean(merged$coverage),

      control_mean_imrs_z = ctrl_mean,
      delivery_mean_imrs_z = del_mean,
      delta_mean_imrs_z = delta_mean,

      control_median_imrs_z = ctrl_median,
      delivery_median_imrs_z = del_median,
      delta_median_imrs_z = delta_median,

      control_sd_imrs_z = safe_sd(ctrl$imrs_z),
      delivery_sd_imrs_z = safe_sd(del$imrs_z),

      cohens_d_imrs_z = d_val,
      t_test_p_imrs_z = t_p,

      auc_imrs_z = auc_val,
      auc_note = auc_note
    )

    sample_level_rows[[length(sample_level_rows) + 1]] <- merged %>%
      mutate(
        gse_id = ds,
        dataset_id = sample_id,
        eval_time_h = eval_time_h
      ) %>%
      transmute(
        gse_id = gse_id,
        dataset_id = dataset_id,
        split_id = split_id,
        contrast_label = contrast_label,
        control_label = control_label,
        tissue = tissue,
        time_h = eval_time_h,
        sample_id = sample_id,
        condition_simple = condition_simple,
        imrs_raw = imrs_raw,
        imrs_z = imrs_z,
        coverage = coverage
      )
  }
}

eval_tbl <- bind_rows(eval_rows)
sample_tbl <- bind_rows(sample_level_rows)

if (nrow(eval_tbl) == 0) {
  stop("No evaluation rows were generated.")
}

summary_tbl <- eval_tbl %>%
  filter(pass) %>%
  group_by(gse_id, tissue) %>%
  summarise(
    n_contrasts = n(),
    mean_delta_mean_imrs_z = safe_mean(delta_mean_imrs_z),
    median_delta_mean_imrs_z = safe_median(delta_mean_imrs_z),
    mean_auc_imrs_z = safe_mean(auc_imrs_z),
    .groups = "drop"
  )

# -------------------------
# WRITE TABLES
# -------------------------
out_eval    <- file.path(eval_dir, "step09_split_eval.tsv")
out_summary <- file.path(eval_dir, "step09_split_summary.tsv")
out_samples <- file.path(eval_dir, "step09_split_sample_level.tsv")

write_tsv(eval_tbl, out_eval)
write_tsv(summary_tbl, out_summary)
write_tsv(sample_tbl, out_samples)

# -------------------------
# PLOTTING
# -------------------------
plot_df <- eval_tbl %>% filter(pass)
sample_plot_df <- sample_tbl %>% filter(is.finite(imrs_z))

if (nrow(plot_df) > 0) {

  # 1. Delta distribution
  p1 <- ggplot(plot_df, aes(x = delta_mean_imrs_z)) +
    geom_histogram(bins = 40) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      title = "IMRS calibration delta distribution",
      x = "Delta mean IMRS_z (DELIVERY - CONTROL)",
      y = "Number of contrasts"
    )
  save_plot_both(p1, file.path(fig_dir, "01_delta_distribution"), 8, 6)

  # 2. Delta by tissue
  p2 <- ggplot(plot_df, aes(x = tissue, y = delta_mean_imrs_z)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      title = "IMRS calibration delta by tissue",
      x = "Tissue",
      y = "Delta mean IMRS_z"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot_both(p2, file.path(fig_dir, "02_delta_by_tissue"), 9, 6)

  # 3. AUC distribution
  auc_plot_df <- plot_df %>% filter(is.finite(auc_imrs_z))
  if (nrow(auc_plot_df) > 0) {
    p3 <- ggplot(auc_plot_df, aes(x = auc_imrs_z)) +
      geom_histogram(bins = 30) +
      geom_vline(xintercept = 0.5, linetype = "dashed") +
      theme_bw() +
      labs(
        title = "IMRS calibration AUC distribution",
        x = "AUC",
        y = "Number of contrasts"
      )
    save_plot_both(p3, file.path(fig_dir, "03_auc_distribution"), 8, 6)
  }

  # 4. Delta vs time
  time_plot_df <- plot_df %>%
    filter(is.finite(time_h), is.finite(delta_mean_imrs_z))

  if (nrow(time_plot_df) > 0) {
    p4 <- ggplot(time_plot_df, aes(x = time_h, y = delta_mean_imrs_z)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw() +
      labs(
        title = "IMRS response over time",
        x = "Hours post delivery",
        y = "Delta mean IMRS_z"
      )

    # add smooth only when enough x diversity exists
    if (length(unique(time_plot_df$time_h)) >= 3) {
      p4 <- p4 + geom_smooth(method = "loess", se = TRUE)
    }

    save_plot_both(p4, file.path(fig_dir, "04_delta_vs_time"), 8, 6)
  }

  # 5. Volcano-style plot
  volcano_df <- plot_df %>%
    mutate(neg_log10_p = -log10(pmax(t_test_p_imrs_z, 1e-300))) %>%
    filter(is.finite(delta_mean_imrs_z), is.finite(neg_log10_p))

  if (nrow(volcano_df) > 0) {
    p5 <- ggplot(volcano_df, aes(x = delta_mean_imrs_z, y = neg_log10_p)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_bw() +
      labs(
        title = "Calibration volcano-style plot",
        x = "Delta mean IMRS_z",
        y = "-log10(t-test p)"
      )
    save_plot_both(p5, file.path(fig_dir, "05_volcano_delta_vs_p"), 8, 6)
  }

  # 6. IMRS_z by condition
  if (nrow(sample_plot_df) > 0) {
    p6 <- ggplot(sample_plot_df, aes(x = condition_simple, y = imrs_z)) +
      geom_boxplot() +
      theme_bw() +
      labs(
        title = "IMRS_z by condition across passing contrasts",
        x = "Condition",
        y = "IMRS_z"
      )
    save_plot_both(p6, file.path(fig_dir, "06_imrsz_by_condition"), 7, 6)
  }

  # 7. Top contrasts barplot
  top_df <- plot_df %>%
    arrange(desc(delta_mean_imrs_z)) %>%
    head(20)

  if (nrow(top_df) > 0) {
    p7 <- ggplot(
      top_df,
      aes(x = reorder(split_id, delta_mean_imrs_z), y = delta_mean_imrs_z)
    ) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      labs(
        title = "Top calibration contrasts by delta IMRS_z",
        x = "Split contrast",
        y = "Delta mean IMRS_z"
      )
    save_plot_both(p7, file.path(fig_dir, "07_top_contrasts_delta"), 10, 8)
  }

  # 8. Delta by GSE
  p8 <- ggplot(plot_df, aes(x = gse_id, y = delta_mean_imrs_z)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      title = "Delta mean IMRS_z by GSE",
      x = "GSE",
      y = "Delta mean IMRS_z"
    )
  save_plot_both(p8, file.path(fig_dir, "08_delta_by_gse"), 8, 6)

  # 9. Effect size vs delta
  effect_df <- plot_df %>%
    filter(is.finite(cohens_d_imrs_z), is.finite(delta_mean_imrs_z))

  if (nrow(effect_df) > 0) {
    p9 <- ggplot(effect_df, aes(x = delta_mean_imrs_z, y = cohens_d_imrs_z)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_bw() +
      labs(
        title = "Effect size versus delta",
        x = "Delta mean IMRS_z",
        y = "Cohen's d"
      )
    save_plot_both(p9, file.path(fig_dir, "09_effect_size_vs_delta"), 8, 6)
  }
}

# -------------------------
# FINISH
# -------------------------
message("\nDONE: Step 09 split-based evaluation finished.")
message("Detailed table: ", out_eval)
message("Summary table:  ", out_summary)
message("Sample table:   ", out_samples)
message("Figures dir:    ", fig_dir)

message("\nTop passing contrasts:")
print(
  eval_tbl %>%
    filter(pass) %>%
    arrange(desc(delta_mean_imrs_z)) %>%
    select(gse_id, split_id, tissue, time_h, n_controls, n_delivery,
           delta_mean_imrs_z, auc_imrs_z) %>%
    head(20)
)