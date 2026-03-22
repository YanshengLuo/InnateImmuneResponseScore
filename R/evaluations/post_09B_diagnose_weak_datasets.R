#!/usr/bin/env Rscript

# ============================================================
# IMRS Step 10 — Diagnose weak / near-zero datasets
#
# Expected inputs:
#   <project_root>/05_score/failure_diagnosis/low_imrs_datasets.tsv
#   <project_root>/05_score/transfer/scores/<ID>__imrs_scores.tsv
#   <project_root>/05_score/transfer/qc/<ID>__qc_summary.tsv
#   <project_root>/05_score/transfer/qc/<ID>__top_contributors.tsv
#   <project_root>/00_metadata/verified_metadata/scoring/<DATASET>/<DATASET>_design.tsv
#
# Optional:
#   DE result table if present somewhere under 04_de/<DATASET>/ or results/<DATASET>/de/
#
# Outputs:
#   <project_root>/05_score/failure_diagnosis/diagnosis/
#       diagnosis_table.tsv
#       diagnosis_report.txt
#       extracted/<ID>/... copied input files
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

failure_root <- file.path(project_root, "05_score", "failure_diagnosis")
low_file     <- file.path(failure_root, "low_imrs_datasets.tsv")
diag_dir     <- file.path(failure_root, "diagnosis")
extract_dir  <- file.path(diag_dir, "extracted")

dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

# thresholds
high_variance_threshold       <- 3.0
control_spread_threshold      <- 1.5
top1_dominance_threshold      <- 0.50
top5_dominance_threshold      <- 0.80
direction_consistency_low     <- 0.60
gene_activation_low_threshold <- 0.30
reversed_score_threshold      <- 0.00

# -------------------------
# HELPERS
# -------------------------
read_table_robust <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) return(NULL)

  x <- tryCatch(read_tsv(path, show_col_types = FALSE, progress = FALSE),
                error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 1) return(x)

  x <- tryCatch(read_csv(path, show_col_types = FALSE, progress = FALSE),
                error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 1) return(x)

  x <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
                error = function(e) NULL)
  if (!is.null(x) && ncol(x) >= 1) return(x)

  NULL
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

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

fmt_chr <- function(x) {
  ifelse(is.na(x) | x == "", "NA", as.character(x))
}

copy_if_exists <- function(src, dst) {
  if (!is.na(src) && src != "" && file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    return(TRUE)
  }
  FALSE
}

normalize_condition_labels <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("CTRL", "UNTREATED", "PBS", "VEHICLE", "NAIVE")] <- "CONTROL"
  x[x %in% c("TREATED", "TX")] <- "DELIVERY"
  x
}

infer_condition_column <- function(df) {
  for (cand in c("condition_simple", "condition", "group", "treatment_group")) {
    if (cand %in% names(df)) return(cand)
  }
  NA_character_
}

infer_contrib_abs_col <- function(df) {
  if ("abs_contribution" %in% names(df)) return("abs_contribution")
  if ("contribution_abs" %in% names(df)) return("contribution_abs")
  NA_character_
}

infer_contrib_col <- function(df) {
  if ("contribution" %in% names(df)) return("contribution")
  if ("weighted_z" %in% names(df)) return("weighted_z")
  NA_character_
}

infer_rank_col <- function(df) {
  if ("rank" %in% names(df)) return("rank")
  if ("contrib_rank" %in% names(df)) return("contrib_rank")
  NA_character_
}

infer_gene_col <- function(df) {
  for (cand in c("gene_id", "gene", "symbol")) {
    if (cand %in% names(df)) return(cand)
  }
  NA_character_
}

infer_log2fc_col <- function(df) {
  for (cand in c("log2FC", "log2fc", "logFC", "beta", "effect", "estimate")) {
    if (cand %in% names(df)) return(cand)
  }
  NA_character_
}

infer_fdr_col <- function(df) {
  for (cand in c("FDR", "padj", "adj_p_val", "qvalue")) {
    if (cand %in% names(df)) return(cand)
  }
  NA_character_
}

find_de_file <- function(project_root, dataset_id) {
  candidates <- character(0)

  roots_to_check <- c(
    file.path(project_root, "04_de", dataset_id),
    file.path(project_root, "results", dataset_id, "de"),
    file.path(project_root, "05_score", "transfer", "de", dataset_id)
  )

  for (rt in roots_to_check) {
    if (dir.exists(rt)) {
      hits <- list.files(
        rt,
        pattern = "(delivery_vs_control|de).*\\.(tsv|csv)$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
      candidates <- c(candidates, hits)
    }
  }

  if (length(candidates) == 0) return(NA_character_)

  preferred <- candidates[grepl("delivery_vs_control", basename(candidates), ignore.case = TRUE)]
  if (length(preferred) > 0) return(preferred[1])

  candidates[1]
}

classify_diagnosis <- function(qc_pass, fail_reason,
                               delivery_mean, delivery_sd, control_sd,
                               top1_median_frac, top5_median_frac,
                               direction_consistency_median,
                               gene_activation_fraction) {
  diagnosis <- "needs_manual_review"
  diagnosis_reason <- character(0)

  if (!is.na(qc_pass) && !qc_pass) {
    diagnosis <- "qc_or_input_issue"
    diagnosis_reason <- c(diagnosis_reason, paste0("QC failed: ", fmt_chr(fail_reason)))
  } else {
    if (!is.na(delivery_mean) && delivery_mean < reversed_score_threshold) {
      diagnosis <- "possible_design_issue"
      diagnosis_reason <- c(diagnosis_reason, "Delivery mean IMRS_z is negative")
    }

    if (!is.na(control_sd) && control_sd > control_spread_threshold) {
      diagnosis <- "possible_design_issue"
      diagnosis_reason <- c(diagnosis_reason, paste0("Control spread high (SD=", fmt_num(control_sd, 2), ")"))
    }

    if (!is.na(delivery_sd) && delivery_sd > high_variance_threshold) {
      diagnosis <- "heterogeneous_response"
      diagnosis_reason <- c(diagnosis_reason, paste0("Delivery variance high (SD=", fmt_num(delivery_sd, 2), ")"))
    }

    if (!is.na(top1_median_frac) && top1_median_frac >= top1_dominance_threshold) {
      diagnosis <- "score_degeneracy"
      diagnosis_reason <- c(diagnosis_reason, paste0("Top1 dominance high (", fmt_num(top1_median_frac), ")"))
    }

    if (!is.na(top5_median_frac) && top5_median_frac >= top5_dominance_threshold) {
      diagnosis <- "score_degeneracy"
      diagnosis_reason <- c(diagnosis_reason, paste0("Top5 dominance high (", fmt_num(top5_median_frac), ")"))
    }

    if (!is.na(direction_consistency_median) && direction_consistency_median < direction_consistency_low) {
      if (diagnosis %in% c("needs_manual_review", "true_weak_biology")) {
        diagnosis <- "heterogeneous_response"
      }
      diagnosis_reason <- c(diagnosis_reason, paste0("Contributor sign consistency low (", fmt_num(direction_consistency_median), ")"))
    }

    if (!is.na(gene_activation_fraction)) {
      if (gene_activation_fraction < gene_activation_low_threshold &&
          (is.na(delivery_sd) || delivery_sd <= high_variance_threshold) &&
          (is.na(top1_median_frac) || top1_median_frac < top1_dominance_threshold)) {
        diagnosis <- "true_weak_biology"
        diagnosis_reason <- c(diagnosis_reason, paste0("Positive-gene fraction low (", fmt_num(gene_activation_fraction), ")"))
      }
    } else {
      if (!is.na(delivery_mean) &&
          delivery_mean >= 0 &&
          delivery_mean <= 10 &&
          (is.na(delivery_sd) || delivery_sd <= high_variance_threshold) &&
          (is.na(top1_median_frac) || top1_median_frac < top1_dominance_threshold) &&
          diagnosis == "needs_manual_review") {
        diagnosis <- "true_weak_biology"
        diagnosis_reason <- c(diagnosis_reason, "Weak positive score without strong variance or dominance signal")
      }
    }
  }

  if (length(diagnosis_reason) == 0) {
    diagnosis_reason <- "No single dominant issue from automated diagnostics"
  } else {
    diagnosis_reason <- paste(unique(diagnosis_reason), collapse = "; ")
  }

  list(
    diagnosis = diagnosis,
    diagnosis_reason = diagnosis_reason
  )
}

# -------------------------
# CORE FUNCTION
# -------------------------
diagnose_dataset <- function(row, project_root, extract_dir) {
  dataset_id <- if ("dataset_id" %in% names(row)) as.character(row[["dataset_id"]][1]) else NA_character_
  id         <- if ("id" %in% names(row)) as.character(row[["id"]][1]) else NA_character_
  score_file <- if ("score_file" %in% names(row)) as.character(row[["score_file"]][1]) else NA_character_
  qc_file    <- if ("qc_file" %in% names(row)) as.character(row[["qc_file"]][1]) else NA_character_

  design_file <- if ("design_file" %in% names(row)) {
    as.character(row[["design_file"]][1])
  } else {
    file.path(project_root, "00_metadata", "verified_metadata", "scoring", dataset_id, paste0(dataset_id, "_design.tsv"))
  }

  class_label <- if ("class" %in% names(row)) as.character(row[["class"]][1]) else NA_character_

  if (is.na(dataset_id) || dataset_id == "") stop("diagnose_dataset(): dataset_id is missing")
  if (is.na(id) || id == "") stop("diagnose_dataset(): id is missing")
  if (is.na(score_file) || score_file == "") stop("diagnose_dataset(): score_file is missing")

  top_file <- file.path(project_root, "05_score", "transfer", "qc", paste0(id, "__top_contributors.tsv"))
  de_file  <- find_de_file(project_root, dataset_id)

  ds_extract_dir <- file.path(extract_dir, make.names(id))
  dir.create(ds_extract_dir, recursive = TRUE, showWarnings = FALSE)

  copy_if_exists(score_file,  file.path(ds_extract_dir, basename(score_file)))
  copy_if_exists(qc_file,     file.path(ds_extract_dir, basename(qc_file)))
  copy_if_exists(top_file,    file.path(ds_extract_dir, basename(top_file)))
  copy_if_exists(design_file, file.path(ds_extract_dir, basename(design_file)))
  copy_if_exists(de_file,     file.path(ds_extract_dir, basename(de_file)))

  scores <- read_table_robust(score_file)
  qc     <- read_table_robust(qc_file)
  top    <- read_table_robust(top_file)
  design <- read_table_robust(design_file)
  de     <- read_table_robust(de_file)

  # score-level metrics
  control_mean  <- NA_real_
  control_sd    <- NA_real_
  delivery_mean <- NA_real_
  delivery_sd   <- NA_real_
  delta_imrs    <- NA_real_
  n_control     <- NA_integer_
  n_delivery    <- NA_integer_

  if (!is.null(scores) && all(c("condition_simple", "imrs_z") %in% names(scores))) {
    scores <- scores %>%
      mutate(condition_simple = normalize_condition_labels(condition_simple))

    ctrl <- scores %>% filter(condition_simple == "CONTROL")
    del  <- scores %>% filter(condition_simple == "DELIVERY")

    control_mean  <- safe_mean(ctrl$imrs_z)
    control_sd    <- safe_sd(ctrl$imrs_z)
    delivery_mean <- safe_mean(del$imrs_z)
    delivery_sd   <- safe_sd(del$imrs_z)
    delta_imrs    <- ifelse(is.na(delivery_mean) | is.na(control_mean), NA_real_, delivery_mean - control_mean)
    n_control     <- nrow(ctrl)
    n_delivery    <- nrow(del)

  } else if (!is.null(design) && !is.null(scores) && "sample_id" %in% names(scores) && "imrs_z" %in% names(scores)) {
    cond_col <- infer_condition_column(design)
    if ("sample_id" %in% names(design) && !is.na(cond_col)) {
      tmp <- scores %>%
        left_join(
          design %>%
            transmute(
              sample_id = as.character(sample_id),
              condition_simple = normalize_condition_labels(.data[[cond_col]])
            ),
          by = "sample_id"
        )

      ctrl <- tmp %>% filter(condition_simple == "CONTROL")
      del  <- tmp %>% filter(condition_simple == "DELIVERY")

      control_mean  <- safe_mean(ctrl$imrs_z)
      control_sd    <- safe_sd(ctrl$imrs_z)
      delivery_mean <- safe_mean(del$imrs_z)
      delivery_sd   <- safe_sd(del$imrs_z)
      delta_imrs    <- ifelse(is.na(delivery_mean) | is.na(control_mean), NA_real_, delivery_mean - control_mean)
      n_control     <- nrow(ctrl)
      n_delivery    <- nrow(del)
    }
  }

  # QC metrics
  qc_pass <- NA
  fail_reason <- NA_character_
  coverage <- NA_real_
  n_genes_sd_floor <- NA_real_

  if (!is.null(qc) && nrow(qc) >= 1) {
    if ("pass" %in% names(qc)) qc_pass <- as.logical(qc$pass[[1]])
    if ("fail_reason" %in% names(qc)) fail_reason <- as.character(qc$fail_reason[[1]])
    if ("coverage" %in% names(qc)) coverage <- suppressWarnings(as.numeric(qc$coverage[[1]]))
    if ("n_genes_sd_floor" %in% names(qc)) n_genes_sd_floor <- suppressWarnings(as.numeric(qc$n_genes_sd_floor[[1]]))
  }

  # contributor metrics
  top1_median_frac <- NA_real_
  top5_median_frac <- NA_real_
  direction_consistency_median <- NA_real_
  contributor_samples_n <- NA_integer_
  top_gene_repeat_fraction <- NA_real_

  if (!is.null(top) && nrow(top) > 0 && "sample_id" %in% names(top)) {
    abs_col  <- infer_contrib_abs_col(top)
    raw_col  <- infer_contrib_col(top)
    rank_col <- infer_rank_col(top)
    gene_col <- infer_gene_col(top)

    if (is.na(abs_col) && !is.na(raw_col)) {
      top <- top %>% mutate(abs_contribution_tmp = abs(.data[[raw_col]]))
      abs_col <- "abs_contribution_tmp"
    }

    if (!is.na(abs_col)) {
      top_work <- top

      if (!is.na(rank_col)) {
        top_work <- top_work %>%
          mutate(.rank_use = suppressWarnings(as.numeric(.data[[rank_col]])))
      } else {
        top_work <- top_work %>%
          group_by(sample_id) %>%
          arrange(desc(.data[[abs_col]]), .by_group = TRUE) %>%
          mutate(.rank_use = row_number()) %>%
          ungroup()
      }

      per_sample <- top_work %>%
        group_by(sample_id) %>%
        summarise(
          top1_frac = {
            x <- .data[[abs_col]]
            s <- sum(x, na.rm = TRUE)
            ifelse(is.finite(s) && s > 0, sum(x[.rank_use <= 1], na.rm = TRUE) / s, NA_real_)
          },
          top5_frac = {
            x <- .data[[abs_col]]
            s <- sum(x, na.rm = TRUE)
            ifelse(is.finite(s) && s > 0, sum(x[.rank_use <= 5], na.rm = TRUE) / s, NA_real_)
          },
          dir_consistency = if (!is.na(raw_col)) {
            mean(sign(.data[[raw_col]]) == sign(sum(.data[[raw_col]], na.rm = TRUE)), na.rm = TRUE)
          } else {
            NA_real_
          },
          .groups = "drop"
        )

      top1_median_frac <- safe_mean(per_sample$top1_frac)
      top5_median_frac <- safe_mean(per_sample$top5_frac)
      direction_consistency_median <- safe_mean(per_sample$dir_consistency)
      contributor_samples_n <- nrow(per_sample)

      if (!is.na(gene_col)) {
        top_gene_repeat_fraction <- top_work %>%
          group_by(sample_id) %>%
          arrange(.rank_use, .by_group = TRUE) %>%
          summarise(top_gene = as.character(.data[[gene_col]][which.min(.rank_use)]), .groups = "drop") %>%
          count(top_gene, name = "n") %>%
          summarise(frac = max(n) / sum(n)) %>%
          pull(frac)

        if (length(top_gene_repeat_fraction) == 0) top_gene_repeat_fraction <- NA_real_
      }
    }
  }

  # DE metrics
  de_found <- !is.null(de)
  de_file_used <- ifelse(is.na(de_file), NA_character_, de_file)
  gene_activation_fraction <- NA_real_
  gene_activation_fraction_fdr05 <- NA_real_
  median_core_log2fc <- NA_real_
  negative_gene_fraction <- NA_real_
  de_gene_n <- NA_integer_

  if (!is.null(de) && nrow(de) > 0) {
    lfc_col <- infer_log2fc_col(de)
    fdr_col <- infer_fdr_col(de)

    if (!is.na(lfc_col)) {
      de[[lfc_col]] <- suppressWarnings(as.numeric(de[[lfc_col]]))
      de_gene_n <- sum(!is.na(de[[lfc_col]]))
      gene_activation_fraction <- mean(de[[lfc_col]] > 0, na.rm = TRUE)
      negative_gene_fraction   <- mean(de[[lfc_col]] < 0, na.rm = TRUE)
      median_core_log2fc       <- stats::median(de[[lfc_col]], na.rm = TRUE)

      if (!is.na(fdr_col)) {
        de[[fdr_col]] <- suppressWarnings(as.numeric(de[[fdr_col]]))
        gene_activation_fraction_fdr05 <- mean(de[[lfc_col]] > 0 & de[[fdr_col]] <= 0.05, na.rm = TRUE)
      }
    }
  }

  diag_res <- classify_diagnosis(
    qc_pass = qc_pass,
    fail_reason = fail_reason,
    delivery_mean = delivery_mean,
    delivery_sd = delivery_sd,
    control_sd = control_sd,
    top1_median_frac = top1_median_frac,
    top5_median_frac = top5_median_frac,
    direction_consistency_median = direction_consistency_median,
    gene_activation_fraction = gene_activation_fraction
  )

  tibble(
    dataset_id = dataset_id,
    id = id,
    class = class_label,
    diagnosis = diag_res$diagnosis,
    diagnosis_reason = diag_res$diagnosis_reason,

    n_control = n_control,
    n_delivery = n_delivery,
    control_mean_imrs_z = control_mean,
    control_sd_imrs_z = control_sd,
    delivery_mean_imrs_z = delivery_mean,
    delivery_sd_imrs_z = delivery_sd,
    delta_imrs = delta_imrs,

    qc_pass = qc_pass,
    fail_reason = fail_reason,
    coverage = coverage,
    n_genes_sd_floor = n_genes_sd_floor,

    top1_median_fraction = top1_median_frac,
    top5_median_fraction = top5_median_frac,
    direction_consistency_median = direction_consistency_median,
    contributor_samples_n = contributor_samples_n,
    top_gene_repeat_fraction = top_gene_repeat_fraction,

    de_found = de_found,
    de_file = de_file_used,
    de_gene_n = de_gene_n,
    gene_activation_fraction = gene_activation_fraction,
    gene_activation_fraction_fdr05 = gene_activation_fraction_fdr05,
    median_core_log2fc = median_core_log2fc,
    negative_gene_fraction = negative_gene_fraction,

    score_file = score_file,
    qc_file = qc_file,
    top_contributors_file = top_file,
    design_file = design_file
  )
}

# -------------------------
# MAIN
# -------------------------
if (!file.exists(low_file)) {
  stop("Missing low_imrs_datasets.tsv: ", low_file)
}

low_tbl <- read_tsv(low_file, show_col_types = FALSE, progress = FALSE)

required_cols <- c("dataset_id", "id", "score_file")
missing_cols <- setdiff(required_cols, names(low_tbl))
if (length(missing_cols) > 0) {
  stop("low_imrs_datasets.tsv missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (nrow(low_tbl) == 0) {
  diag_tbl <- tibble()
} else {
  diag_tbl <- purrr::map_dfr(seq_len(nrow(low_tbl)), function(i) {
    row <- low_tbl[i, , drop = FALSE]
    diagnose_dataset(row, project_root = project_root, extract_dir = extract_dir)
  })
}

diag_tbl <- diag_tbl %>%
  arrange(
    factor(class, levels = c("reversed", "near_zero", "weak", "strong", "unknown")),
    factor(diagnosis, levels = c(
      "qc_or_input_issue",
      "possible_design_issue",
      "score_degeneracy",
      "heterogeneous_response",
      "true_weak_biology",
      "needs_manual_review"
    )),
    delivery_mean_imrs_z
  )

out_tsv <- file.path(diag_dir, "diagnosis_table.tsv")
write_tsv(diag_tbl, out_tsv)

# -------------------------
# REPORT
# -------------------------
report_file <- file.path(diag_dir, "diagnosis_report.txt")

diag_counts <- if (nrow(diag_tbl) > 0) diag_tbl %>% count(diagnosis, name = "n") else tibble(diagnosis = character(), n = integer())
class_counts <- if (nrow(diag_tbl) > 0) diag_tbl %>% count(class, name = "n") else tibble(class = character(), n = integer())

lines_out <- c(
  "IMRS WEAK-DATASET DIAGNOSIS REPORT",
  "==================================",
  "",
  paste0("Project root: ", project_root),
  paste0("Input low table: ", low_file),
  paste0("Output diagnosis table: ", out_tsv),
  "",
  "Interpretation of diagnosis labels:",
  "  qc_or_input_issue     = failed QC or insufficient input quality",
  "  possible_design_issue = reversed signal or suspicious control spread",
  "  score_degeneracy      = score dominated by one/few genes",
  "  heterogeneous_response = heterogeneous response across samples",
  "  true_weak_biology     = genuinely weak signal without strong noise/dominance flags",
  "  needs_manual_review   = automated metrics inconclusive",
  "",
  paste0("Flagged datasets analyzed: ", nrow(diag_tbl)),
  ""
)

lines_out <- c(lines_out, "Counts by class:")
if (nrow(class_counts) == 0) {
  lines_out <- c(lines_out, "  (none)")
} else {
  for (i in seq_len(nrow(class_counts))) {
    lines_out <- c(lines_out, paste0("  - ", class_counts$class[i], ": ", class_counts$n[i]))
  }
}

lines_out <- c(lines_out, "", "Counts by diagnosis:")
if (nrow(diag_counts) == 0) {
  lines_out <- c(lines_out, "  (none)")
} else {
  for (i in seq_len(nrow(diag_counts))) {
    lines_out <- c(lines_out, paste0("  - ", diag_counts$diagnosis[i], ": ", diag_counts$n[i]))
  }
}

lines_out <- c(lines_out, "", "Per-dataset summary:", "--------------------")

if (nrow(diag_tbl) == 0) {
  lines_out <- c(lines_out, "  No flagged datasets to diagnose.")
} else {
  for (i in seq_len(nrow(diag_tbl))) {
    x <- diag_tbl[i, ]

    lines_out <- c(
      lines_out,
      paste0("Dataset: ", x$dataset_id, " | id=", x$id),
      paste0("  class=", fmt_chr(x$class), " | diagnosis=", fmt_chr(x$diagnosis)),
      paste0("  delivery_mean_imrs_z=", fmt_num(x$delivery_mean_imrs_z),
             " | control_mean_imrs_z=", fmt_num(x$control_mean_imrs_z),
             " | delta_imrs=", fmt_num(x$delta_imrs)),
      paste0("  delivery_sd=", fmt_num(x$delivery_sd_imrs_z),
             " | control_sd=", fmt_num(x$control_sd_imrs_z),
             " | n_control=", fmt_chr(x$n_control),
             " | n_delivery=", fmt_chr(x$n_delivery)),
      paste0("  qc_pass=", fmt_chr(x$qc_pass),
             " | coverage=", fmt_num(x$coverage),
             " | n_genes_sd_floor=", fmt_chr(x$n_genes_sd_floor)),
      paste0("  top1_median_fraction=", fmt_num(x$top1_median_fraction),
             " | top5_median_fraction=", fmt_num(x$top5_median_fraction),
             " | direction_consistency=", fmt_num(x$direction_consistency_median)),
      paste0("  gene_activation_fraction=", fmt_num(x$gene_activation_fraction),
             " | gene_activation_fraction_fdr05=", fmt_num(x$gene_activation_fraction_fdr05),
             " | median_core_log2fc=", fmt_num(x$median_core_log2fc)),
      paste0("  reason=", fmt_chr(x$diagnosis_reason)),
      paste0("  score_file=", fmt_chr(x$score_file)),
      paste0("  qc_file=", fmt_chr(x$qc_file)),
      paste0("  top_contributors_file=", fmt_chr(x$top_contributors_file)),
      paste0("  design_file=", fmt_chr(x$design_file)),
      paste0("  de_file=", fmt_chr(x$de_file)),
      ""
    )
  }
}

writeLines(lines_out, report_file)

message("Saved:")
message("  ", out_tsv)
message("  ", report_file)
message("  Extracted files copied under: ", extract_dir)
message("Done.")