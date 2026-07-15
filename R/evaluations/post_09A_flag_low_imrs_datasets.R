#!/usr/bin/env Rscript

# ============================================================
# IMRS Step 09 — Flag low / weak datasets and generate report
#
# Inputs expected:
#   <project_root>/05_score/transfer/scores/<ID>__imrs_scores.tsv
#   <project_root>/05_score/transfer/qc/<ID>__qc_summary.tsv
#   <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#
# Outputs:
#   <project_root>/05_score/failure_diagnosis/
#       all_dataset_imrs_summary.tsv
#       low_imrs_datasets.tsv
#       low_imrs_review_table.tsv
#       mean_imrs_effect_by_dataset.png
#       low_imrs_report.txt
#
# Classification:
#   reversed  : delivery_mean_imrs_z < 0
#   near_zero : 0 <= delivery_mean_imrs_z <= 0.5
#   weak      : 0.5 < delivery_mean_imrs_z <= 10
#   strong    : delivery_mean_imrs_z > 10
#
# low_imrs_flag includes reversed + near_zero + weak
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

scores_dir  <- file.path(project_root, "05_score", "transfer", "scores")
qc_dir      <- file.path(project_root, "05_score", "transfer", "qc")
design_root <- file.path(project_root, "00_metadata", "verified_metadata", "scoring")
out_dir     <- file.path(project_root, "05_score", "failure_diagnosis")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

near_zero_threshold <- 0.5
weak_threshold      <- 10

# Edit this if you want different anchor labels in the plot/report
anchor_ids <- c("GSE39129", "GSE167521", "GSE279744", "GSE264344")

# -------------------------
# HELPERS
# -------------------------
read_tsv_safe <- function(path) {
  tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
}

extract_id_info <- function(score_path) {
  base <- basename(score_path)
  id <- sub("__imrs_scores\\.tsv$", "", base)
  dataset_id <- str_extract(id, "^GSE[0-9]+")
  tibble(
    id = id,
    dataset_id = dataset_id
  )
}

safe_mean <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

ci95_halfwidth <- function(x) {
  x <- stats::na.omit(x)
  n <- length(x)
  if (n < 2) return(NA_real_)
  s <- sd(x)
  tcrit <- qt(0.975, df = n - 1)
  tcrit * s / sqrt(n)
}

classify_dataset <- function(delivery_mean_imrs_z) {
  if (is.na(delivery_mean_imrs_z)) return("unknown")
  if (delivery_mean_imrs_z < 0) {
    "reversed"
  } else if (delivery_mean_imrs_z <= near_zero_threshold) {
    "near_zero"
  } else if (delivery_mean_imrs_z <= weak_threshold) {
    "weak"
  } else {
    "strong"
  }
}

dataset_type_label <- function(dataset_id) {
  ifelse(dataset_id %in% anchor_ids, "anchor", "external")
}

find_design_file <- function(dataset_id) {
  candidate <- file.path(design_root, dataset_id, paste0(dataset_id, "_design.tsv"))
  if (file.exists(candidate)) candidate else NA_character_
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

fmt_chr <- function(x) {
  ifelse(is.na(x) | x == "", "NA", as.character(x))
}

infer_issue_summary <- function(pass, fail_reason, coverage, n_genes_sd_floor,
                                control_sd_imrs_z, delivery_sd_imrs_z,
                                class_label) {
  reasons <- character(0)

  if (!is.na(pass) && !pass) {
    reasons <- c(reasons, paste0("QC fail: ", fmt_chr(fail_reason)))
  }

  if (!is.na(coverage) && coverage < 0.9) {
    reasons <- c(reasons, paste0("Lower gene coverage (", fmt_num(coverage), ")"))
  }

  if (!is.na(n_genes_sd_floor) && n_genes_sd_floor >= 50) {
    reasons <- c(reasons, paste0("Many genes hit SD floor (n=", n_genes_sd_floor, ")"))
  }

  if (!is.na(delivery_sd_imrs_z) && delivery_sd_imrs_z > 5) {
    reasons <- c(reasons, paste0("High delivery-group variance (SD=", fmt_num(delivery_sd_imrs_z, 2), ")"))
  }

  if (!is.na(control_sd_imrs_z) && control_sd_imrs_z > 2) {
    reasons <- c(reasons, paste0("Control spread larger than expected (SD=", fmt_num(control_sd_imrs_z, 2), ")"))
  }

  if (length(reasons) == 0) {
    if (class_label == "reversed") {
      reasons <- c(reasons, "Reversed direction; check design labels, control grouping, timepoint, or out-of-scope biology")
    } else if (class_label == "near_zero") {
      reasons <- c(reasons, "Near-zero signal; likely weak biology, late timepoint, or reduced transferability")
    } else if (class_label == "weak") {
      reasons <- c(reasons, "Weak positive signal; investigate biological context, variance, and top contributors")
    } else if (class_label == "strong") {
      reasons <- c(reasons, "No obvious technical issue from summary metrics")
    } else {
      reasons <- c(reasons, "Unable to classify from available metrics")
    }
  }

  paste(reasons, collapse = "; ")
}

# -------------------------
# 1) READ SCORE FILES
# -------------------------
score_files <- list.files(scores_dir, pattern = "__imrs_scores\\.tsv$", full.names = TRUE)

if (length(score_files) == 0) {
  stop("No score files found under: ", scores_dir)
}

score_summary <- map_dfr(score_files, function(f) {
  info <- extract_id_info(f)
  sc <- read_tsv_safe(f)

  if (is.null(sc)) {
    return(tibble(
      id = info$id,
      dataset_id = info$dataset_id,
      score_file = f,
      parse_ok = FALSE,
      n_control = NA_integer_,
      n_delivery = NA_integer_,
      control_mean_imrs_z = NA_real_,
      control_sd_imrs_z = NA_real_,
      delivery_mean_imrs_z = NA_real_,
      delivery_sd_imrs_z = NA_real_,
      delta_imrs = NA_real_,
      ci95_halfwidth = NA_real_
    ))
  }

  required_cols <- c("sample_id", "condition_simple", "imrs_z")
  if (!all(required_cols %in% names(sc))) {
    return(tibble(
      id = info$id,
      dataset_id = info$dataset_id,
      score_file = f,
      parse_ok = FALSE,
      n_control = NA_integer_,
      n_delivery = NA_integer_,
      control_mean_imrs_z = NA_real_,
      control_sd_imrs_z = NA_real_,
      delivery_mean_imrs_z = NA_real_,
      delivery_sd_imrs_z = NA_real_,
      delta_imrs = NA_real_,
      ci95_halfwidth = NA_real_
    ))
  }

  sc <- sc %>%
    mutate(condition_simple = toupper(as.character(condition_simple)))

  ctrl <- sc %>% filter(condition_simple == "CONTROL")
  del  <- sc %>% filter(condition_simple == "DELIVERY")

  ctrl_mean <- safe_mean(ctrl$imrs_z)
  ctrl_sd   <- safe_sd(ctrl$imrs_z)
  del_mean  <- safe_mean(del$imrs_z)
  del_sd    <- safe_sd(del$imrs_z)
  delta     <- ifelse(is.na(del_mean) | is.na(ctrl_mean), NA_real_, del_mean - ctrl_mean)
  ci_hw     <- ci95_halfwidth(del$imrs_z)

  tibble(
    id = info$id,
    dataset_id = info$dataset_id,
    score_file = f,
    parse_ok = TRUE,
    n_control = nrow(ctrl),
    n_delivery = nrow(del),
    control_mean_imrs_z = ctrl_mean,
    control_sd_imrs_z = ctrl_sd,
    delivery_mean_imrs_z = del_mean,
    delivery_sd_imrs_z = del_sd,
    delta_imrs = delta,
    ci95_halfwidth = ci_hw
  )
})

# -------------------------
# 2) READ QC FILES
# -------------------------
qc_files <- list.files(qc_dir, pattern = "__qc_summary\\.tsv$", full.names = TRUE)

qc_summary <- map_dfr(qc_files, function(f) {
  q <- read_tsv_safe(f)
  if (is.null(q) || nrow(q) == 0) return(NULL)

  # Rename potentially colliding columns only if they exist
  if ("n_controls" %in% names(q)) {
    q <- q %>% rename(n_controls_qc = n_controls)
  }
  if ("n_delivery" %in% names(q)) {
    q <- q %>% rename(n_delivery_qc = n_delivery)
  }

  q %>%
    mutate(qc_file = f) %>%
    select(any_of(c(
      "id", "dataset_id", "counts_path", "pass", "fail_reason",
      "n_samples", "n_controls_qc", "n_delivery_qc", "n_weights", "n_overlap",
      "coverage", "n_genes_sd_floor", "gene_sd_floor", "score_sd_floor",
      "controls_imrs_z_mean", "controls_imrs_z_sd", "qc_file"
    )))
})

# -------------------------
# 3) MERGE + CLASSIFY
# -------------------------
summary_tbl <- score_summary %>%
  left_join(qc_summary, by = c("id", "dataset_id")) %>%
  mutate(
    n_control_final  = dplyr::coalesce(n_control, n_controls_qc),
    n_delivery_final = dplyr::coalesce(n_delivery, n_delivery_qc),
    dataset_type     = dataset_type_label(dataset_id),
    class            = vapply(delivery_mean_imrs_z, classify_dataset, character(1)),
    low_imrs_flag    = class %in% c("reversed", "near_zero", "weak"),
    design_file      = vapply(dataset_id, find_design_file, character(1)),
    design_exists    = !is.na(design_file),
    issue_summary    = pmap_chr(
      list(pass, fail_reason, coverage, n_genes_sd_floor,
           control_sd_imrs_z, delivery_sd_imrs_z, class),
      infer_issue_summary
    )
  ) %>%
  arrange(
    factor(class, levels = c("reversed", "near_zero", "weak", "strong", "unknown")),
    delivery_mean_imrs_z
  )

low_tbl <- summary_tbl %>%
  filter(low_imrs_flag)

review_tbl <- low_tbl %>%
  select(
    dataset_id, id, dataset_type, class,
    n_control_final, n_delivery_final,
    control_mean_imrs_z, control_sd_imrs_z,
    delivery_mean_imrs_z, delivery_sd_imrs_z, delta_imrs, ci95_halfwidth,
    pass, fail_reason, coverage, n_overlap, n_weights, n_genes_sd_floor,
    counts_path, score_file, qc_file, design_file, design_exists, issue_summary
  )

# -------------------------
# 4) SAVE TABLES
# -------------------------
write_tsv(summary_tbl, file.path(out_dir, "all_dataset_imrs_summary.tsv"))
write_tsv(low_tbl, file.path(out_dir, "low_imrs_datasets.tsv"))
write_tsv(review_tbl, file.path(out_dir, "low_imrs_review_table.tsv"))

# -------------------------
# 5) PLOT
# -------------------------
plot_df <- summary_tbl %>%
  filter(!is.na(delivery_mean_imrs_z)) %>%
  mutate(
    dataset_label = dataset_id,
    xmin = delivery_mean_imrs_z - ci95_halfwidth,
    xmax = delivery_mean_imrs_z + ci95_halfwidth
  ) %>%
  arrange(delivery_mean_imrs_z) %>%
  mutate(dataset_label = factor(dataset_label, levels = rev(unique(dataset_label))))

p <- ggplot(plot_df, aes(x = delivery_mean_imrs_z, y = dataset_label, color = dataset_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 0.15, na.rm = TRUE) +
  geom_point(size = 2.4) +
  labs(
    title = "Mean IMRS effect size by dataset",
    x = "Mean Δ IMRS_z with 95% CI",
    y = "Dataset",
    color = "Dataset type"
  ) +
  theme_bw(base_size = 11)

ggsave(
  filename = file.path(out_dir, "mean_imrs_effect_by_dataset.png"),
  plot = p,
  width = 8.5, height = 6, dpi = 300
)

# -------------------------
# 6) TEXT REPORT
# -------------------------
report_path <- file.path(out_dir, "low_imrs_report.txt")

n_total <- nrow(summary_tbl)
n_low   <- nrow(low_tbl)

class_counts <- summary_tbl %>%
  count(class, name = "n") %>%
  arrange(match(class, c("reversed", "near_zero", "weak", "strong", "unknown")))

anchor_summary <- summary_tbl %>%
  filter(dataset_type == "anchor") %>%
  select(dataset_id, class, delivery_mean_imrs_z, delta_imrs)

external_summary <- summary_tbl %>%
  filter(dataset_type == "external") %>%
  select(dataset_id, class, delivery_mean_imrs_z, delta_imrs)

lines_out <- c(
  "IMRS FAILURE / WEAK-SIGNAL REPORT",
  "================================",
  "",
  paste0("Project root: ", project_root),
  paste0("Scores dir:   ", scores_dir),
  paste0("QC dir:       ", qc_dir),
  paste0("Design root:  ", design_root),
  "",
  "Classification thresholds:",
  paste0("  reversed  : delivery_mean_imrs_z < 0"),
  paste0("  near_zero : 0 <= delivery_mean_imrs_z <= ", near_zero_threshold),
  paste0("  weak      : ", near_zero_threshold, " < delivery_mean_imrs_z <= ", weak_threshold),
  paste0("  strong    : delivery_mean_imrs_z > ", weak_threshold),
  "",
  paste0("Datasets analyzed: ", n_total),
  paste0("Flagged low/weak datasets: ", n_low),
  "",
  "Class counts:"
)

for (i in seq_len(nrow(class_counts))) {
  lines_out <- c(lines_out, paste0("  - ", class_counts$class[i], ": ", class_counts$n[i]))
}

lines_out <- c(lines_out, "", "Anchor datasets:", "----------------")
if (nrow(anchor_summary) == 0) {
  lines_out <- c(lines_out, "  (none found)")
} else {
  for (i in seq_len(nrow(anchor_summary))) {
    lines_out <- c(
      lines_out,
      paste0(
        "  - ", anchor_summary$dataset_id[i],
        " | class=", anchor_summary$class[i],
        " | delivery_mean_imrs_z=", fmt_num(anchor_summary$delivery_mean_imrs_z[i]),
        " | delta_imrs=", fmt_num(anchor_summary$delta_imrs[i])
      )
    )
  }
}

lines_out <- c(lines_out, "", "External datasets:", "------------------")
if (nrow(external_summary) == 0) {
  lines_out <- c(lines_out, "  (none found)")
} else {
  for (i in seq_len(nrow(external_summary))) {
    lines_out <- c(
      lines_out,
      paste0(
        "  - ", external_summary$dataset_id[i],
        " | class=", external_summary$class[i],
        " | delivery_mean_imrs_z=", fmt_num(external_summary$delivery_mean_imrs_z[i]),
        " | delta_imrs=", fmt_num(external_summary$delta_imrs[i])
      )
    )
  }
}

lines_out <- c(lines_out, "", "Flagged datasets:", "-----------------")
if (nrow(review_tbl) == 0) {
  lines_out <- c(lines_out, "  None flagged under current thresholds.")
} else {
  for (i in seq_len(nrow(review_tbl))) {
    row <- review_tbl[i, ]
    lines_out <- c(
      lines_out,
      paste0("  - ", row$dataset_id,
             " | id=", row$id,
             " | class=", row$class,
             " | type=", row$dataset_type),
      paste0("      delivery_mean_imrs_z=", fmt_num(row$delivery_mean_imrs_z),
             " | control_mean_imrs_z=", fmt_num(row$control_mean_imrs_z),
             " | delta_imrs=", fmt_num(row$delta_imrs)),
      paste0("      n_control=", fmt_chr(row$n_control_final),
             " | n_delivery=", fmt_chr(row$n_delivery_final),
             " | coverage=", fmt_num(row$coverage)),
      paste0("      QC pass=", fmt_chr(row$pass),
             " | fail_reason=", fmt_chr(row$fail_reason)),
      paste0("      n_genes_sd_floor=", fmt_chr(row$n_genes_sd_floor)),
      paste0("      issue_summary=", fmt_chr(row$issue_summary)),
      paste0("      design_file=", fmt_chr(row$design_file)),
      ""
    )
  }
}

lines_out <- c(
  lines_out,
  "Interpretation notes:",
  "---------------------",
  "1) These are dataset-internal control-centered scores; weak datasets are not automatically invalid.",
  "2) Review design.tsv, per-sample scores, and top contributors before concluding biology vs metadata problems.",
  "3) Weak and near-zero datasets are boundary-condition datasets for transferability, not just failures."
)

writeLines(lines_out, con = report_path)

# -------------------------
# 7) CONSOLE SUMMARY
# -------------------------
message("Saved:")
message("  ", file.path(out_dir, "all_dataset_imrs_summary.tsv"))
message("  ", file.path(out_dir, "low_imrs_datasets.tsv"))
message("  ", file.path(out_dir, "low_imrs_review_table.tsv"))
message("  ", file.path(out_dir, "mean_imrs_effect_by_dataset.png"))
message("  ", report_path)
message("Done.")