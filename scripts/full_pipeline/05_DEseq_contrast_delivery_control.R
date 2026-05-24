# ============================================================
# IMRS Step 5 - DESeq2 contrasts: DELIVERY vs CONTROL (WINDOWS)
#
# Part A: Size-factor normalization (from RAW counts + full design)
#   - Compute size factors separately for:
#       (1) ANCHOR phase     : H <= 24
#       (2) CALIBRATION      : size factors computed on ALL samples (no time filtering)
#   - IMPORTANT: Anchor SFs use ALL controls + time-filtered delivery.
#                Calibration SFs use ALL samples (prevents sparse SF overlap across contrasts).
#     Always include ALL CONTROL samples in size-factor estimation.
#     Only time-filter DELIVERY samples.
#   - Hard-validate counts are integers (no rounding)
#
# Part B: DESeq2 contrasts (from RAW counts + split designs)
#   - Filter split-design TSVs by parsing "__H=<hours>__" from filename
#   - Use RAW counts (integer) + precomputed sizeFactors for the PHASE
#   - Write TWO outputs per contrast:
#       1) full DESeq2 table (debug / trace)
#       2) workflow schema: gene, log2FC, SE, stat, pval, FDR
#
# Inputs:
#   RAW counts:
#     <project_root>/03_counts/<DATASET>/featurecounts/validation/gene_counts_clean.tsv
#   FFull design:
#     <project_root>/00_metadata/verified_metadata/<DATASET>_design.tsv
#   Split designs:
#     <project_root>/00_metadata/verified_metadata/splited/<DATASET>_design/*.tsv
#
# Outputs:
#   Normalized (phase-specific):
#     <project_root>/04_de/<DATASET>/normalized/anchor/size_factors.tsv
#     <project_root>/04_de/<DATASET>/normalized/anchor/gene_counts_normalized.tsv
#     <project_root>/04_de/<DATASET>/normalized/calibration/size_factors.tsv
#     <project_root>/04_de/<DATASET>/normalized/calibration/gene_counts_normalized.tsv
#
#   DE results (phase-specific):
#     <project_root>/04_de/comparison/<DATASET>/deseq2_contrasts/anchor/*.tsv
#     <project_root>/04_de/comparison/<DATASET>/deseq2_contrasts/calibration/*.tsv
#
# Tracker output:
#   <project_root>/04_de/comparison/_tracker/step5_phase_tracker.tsv
#
# Usage:
#   Rscript 05_IMRS_Step5_Windows_Phased.R [project_root] [dataset_id_optional]
# Example:
#   Rscript 05_IMRS_Step5_Windows_Phased.R "D:/IMRS_Project" "GSE167521"
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# -------------------------
# CONFIG
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"
dataset_only <- if (length(args) >= 2) args[2] else NA_character_

de_root <- file.path(project_root, "04_de")

split_design_root <- file.path(
  project_root,
  "00_metadata", "verified_metadata", "splited"
)

out_base <- file.path(project_root, "04_de", "comparison")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

min_n_per_group <- 2
gene_min_total  <- 10

# Phase rules (as requested):
#   anchor:      H <= 24
#   calibration: 24 < H < 72
anchor_hours_max      <- 24
calib_hours_lower_exc <- 24
calib_hours_upper_exc <- 73  # your existing convention: <73 means <=72

run_anchor      <- TRUE
run_calibration <- TRUE

message("Project root: ", project_root)
message("Split design root: ", split_design_root)
message("Output base: ", out_base)
message("ANCHOR rule: H <= ", anchor_hours_max)
message("CALIBRATION rule: H > ", calib_hours_lower_exc, " and H < ", calib_hours_upper_exc)

# -------------------------
# TRACKER (metadata + processing audit)
# -------------------------
TRACKER <- new.env(parent = emptyenv())
TRACKER$rows <- list()

`%||%` <- function(a, b) if (!is.null(a)) a else b

track_add <- function(row_list) {
  TRACKER$rows[[length(TRACKER$rows) + 1]] <<- tibble::as_tibble(row_list)
}

track_write <- function(out_base) {
  if (length(TRACKER$rows) == 0) {
    message("Tracker: no rows to write.")
    return(invisible(NULL))
  }
  df <- dplyr::bind_rows(TRACKER$rows)

  out_dir <- file.path(out_base, "_tracker")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(out_dir, "step5_phase_tracker.tsv")
  readr::write_tsv(df, out_file)
  message("Tracker wrote: ", out_file)

  invisible(df)
}

# -------------------------
# HELPERS
# -------------------------
read_table_robust <- function(path) {
  df <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  df <- tryCatch(read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  df <- tryCatch(read_table(path, show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  stop("Failed to read (delimiter/header issue): ", path)
}

make_counts_matrix <- function(df) {

  gene_col <- dplyr::case_when(
    "gene_id" %in% names(df) ~ "gene_id",
    "Geneid"  %in% names(df) ~ "Geneid",
    "GeneID"  %in% names(df) ~ "GeneID",
    "gene"    %in% names(df) ~ "gene",
    TRUE                     ~ names(df)[1]
  )

  names(df)[names(df) == gene_col] <- "gene_id"

  if (ncol(df) < 3) {
    stop("Counts table has <3 columns (need gene_id + >=2 samples). Columns: ",
         paste(names(df), collapse = ", "))
  }

  sample_cols <- setdiff(names(df), "gene_id")
  gene_ids <- as.character(df$gene_id)

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- gene_ids
  colnames(mat) <- trimws(sample_cols)

  if (anyNA(mat)) stop("NA found in counts matrix after coercion. Check file formatting.")
  mat
}

assert_integer_counts <- function(mat, file_label = "counts") {
  if (any(mat < 0, na.rm = TRUE)) {
    stop(file_label, ": negative values detected. Fix upstream.")
  }
  frac <- abs(mat - round(mat))
  if (any(frac > 1e-8, na.rm = TRUE)) {
    stop(
      file_label, ": non-integer counts detected (fractional values). ",
      "DESeq2 requires raw integer counts. Do NOT use TPM/FPKM/normalized matrices here."
    )
  }
}

is_real_split_design <- function(path) {
  bn <- tolower(basename(path))
  if (!str_detect(bn, "\\.tsv$")) return(FALSE)
  if (str_detect(bn, "contrast_index|group_map")) return(FALSE)
  str_detect(bn, "design__t=")
}

# Parse hours from split design filename like "__H=3__"
parse_H_from_name <- function(path) {
  bn <- basename(path)
  m <- stringr::str_match(bn, "__H=([^_]+)__")
  if (is.na(m[1,2])) return(NA_real_)
  suppressWarnings(as.numeric(m[1,2]))
}

phase_from_H <- function(H,
                         anchor_max = 24,
                         calib_lower_exc = 24,
                         calib_upper_exc = 72) {
  if (is.na(H)) return(NA_character_)
  if (H <= anchor_max) return("anchor")
  if (H > calib_lower_exc && H < calib_upper_exc) return("calibration")
  return("other")
}

# Try to find time hours in FULL design (for SF estimation); returns numeric hours or all-NA
get_time_hours <- function(design_df) {
  nm <- names(design_df)
  nm_low <- tolower(nm)

  cand <- c(
    "time_hr", "time_hrs", "time_hours", "hours", "hour",
    "timepoint_hr", "timepoint_hrs", "timepoint_hours",
    "t_hr", "t_hrs"
  )
  idx <- match(cand, nm_low)
  idx <- idx[!is.na(idx)][1]
  if (!is.na(idx)) {
    col <- nm[idx]
    return(suppressWarnings(as.numeric(design_df[[col]])))
  }

  cand2 <- c("time", "timepoint", "treatment_time", "post_injection_time", "post_infection_time")
  idx2 <- match(cand2, nm_low)
  idx2 <- idx2[!is.na(idx2)][1]
  if (!is.na(idx2)) {
    col <- nm[idx2]
    raw <- as.character(design_df[[col]])
    return(suppressWarnings(as.numeric(str_extract(raw, "[0-9]+\\.?[0-9]*"))))
  }

  rep(NA_real_, nrow(design_df))
}

# Anchors: prefer explicit anchor flag if present, else time <= cutoff if available.
anchor_filter_full_design <- function(design_df, cutoff_hours = 24) {
  nm_low <- tolower(names(design_df))

  if ("anchor" %in% nm_low) {
    col <- names(design_df)[which(nm_low == "anchor")[1]]
    v <- design_df[[col]]
    return(tolower(as.character(v)) %in% c("1","true","t","yes","y","anchor"))
  }
  if ("anchor_flag" %in% nm_low) {
    col <- names(design_df)[which(nm_low == "anchor_flag")[1]]
    v <- design_df[[col]]
    return(tolower(as.character(v)) %in% c("1","true","t","yes","y"))
  }
  if ("phase" %in% nm_low) {
    col <- names(design_df)[which(nm_low == "phase")[1]]
    v <- tolower(as.character(design_df[[col]]))
    keep <- v %in% c("anchor","anchors")
    if (any(keep)) return(keep)
  }

  th <- get_time_hours(design_df)
  if (all(is.na(th))) return(rep(TRUE, nrow(design_df)))
  keep <- !is.na(th) & th <= cutoff_hours
  if (!any(keep)) return(rep(TRUE, nrow(design_df)))
  keep
}

# Try to find CONTROL/DELIVERY labels in FULL design (for SF estimation)
get_condition_simple_full_design <- function(design_df) {
  nm_low <- tolower(names(design_df))

  cand <- c("condition_simple", "condition", "group", "treatment", "arm")
  idx <- match(cand, nm_low)
  idx <- idx[!is.na(idx)][1]
  if (is.na(idx)) return(rep(NA_character_, nrow(design_df)))

  col <- names(design_df)[idx]
  v <- tolower(trimws(as.character(design_df[[col]])))

  out <- rep(NA_character_, length(v))
  out[str_detect(v, "control|baseline|vehicle|mock|untreat|naive|saline|pbs")] <- "CONTROL"
  out[str_detect(v, "delivery|treat|treated|infect|inj|ad5|ad26|aav|lenti|dose|vax|vacc")] <- "DELIVERY"
  out
}

# Metadata tracker: count split designs by phase (mouse-only assumption is handled upstream in how you call datasets)
count_split_designs_by_phase <- function(dataset_id,
                                        anchor_max = 24,
                                        calib_lower_exc = 24,
                                        calib_upper_exc = 72) {
  design_dir <- file.path(split_design_root, paste0(dataset_id, "_design"))
  if (!dir.exists(design_dir)) {
    return(list(total = 0, anchor = 0, calibration = 0, other = 0, missingH = 0))
  }

  design_files <- list.files(design_dir, pattern = "\\.tsv$", full.names = TRUE)
  design_files <- design_files[sapply(design_files, is_real_split_design)]
  if (length(design_files) == 0) {
    return(list(total = 0, anchor = 0, calibration = 0, other = 0, missingH = 0))
  }

  Hs <- vapply(design_files, parse_H_from_name, numeric(1))
  phs <- vapply(
    Hs, phase_from_H, character(1),
    anchor_max = anchor_max,
    calib_lower_exc = calib_lower_exc,
    calib_upper_exc = calib_upper_exc
  )

  list(
    total = length(design_files),
    anchor = sum(phs == "anchor", na.rm = TRUE),
    calibration = sum(phs == "calibration", na.rm = TRUE),
    other = sum(phs == "other", na.rm = TRUE),
    missingH = sum(is.na(Hs))
  )
}

# -------------------------
# STEP 5A: compute phase-specific size factors + normalized counts
# -------------------------
compute_size_factors <- function(dataset_id, raw_mat, full_design_path,
                                 phase = c("anchor","calibration"),
                                 anchor_max = 24,
                                 calib_lower_exc = 24,
                                 calib_upper_exc = 72) {
  phase <- match.arg(phase)

  if (!file.exists(full_design_path)) stop("Full design file not found: ", full_design_path)

  design <- read_table_robust(full_design_path)
  if (!("sample_id" %in% colnames(design))) stop("Full design must contain 'sample_id': ", full_design_path)

  design <- design %>%
    mutate(sample_id = trimws(as.character(sample_id))) %>%
    filter(sample_id %in% colnames(raw_mat))

  if (nrow(design) < 4) stop("Too few samples after aligning full design to counts for ", dataset_id)

  # ---- KEY CHANGE ----
  # For calibration: compute size factors on ALL samples in the dataset
  if (phase == "calibration") {
    design_sf <- design
    message("Size factors computed on phase=calibration for ", dataset_id, " using ALL samples (n=", nrow(design_sf), ").")
  } else {
    # ---- ANCHOR behavior unchanged (ALL controls + time-filtered delivery if possible) ----
    th <- get_time_hours(design)
    has_time <- !all(is.na(th))

    cond_full <- get_condition_simple_full_design(design)
    has_cond  <- !all(is.na(cond_full))

    # If we can't reliably identify condition in FULL design, fall back to ALL samples
    if (!has_cond) {
      message("NOTE: Could not identify CONTROL/DELIVERY in full design for ", dataset_id,
              "; computing size factors on ALL samples for phase=anchor")
      design_sf <- design
    } else {
      is_control  <- cond_full == "CONTROL"
      is_delivery <- cond_full == "DELIVERY"

      if (has_time) {
        keep_delivery <- !is.na(th) & th <= anchor_max
      } else {
        # No time column: fall back to anchor_flag/phase heuristic (DELIVERY only)
        is_anchor <- anchor_filter_full_design(design, cutoff_hours = anchor_max)
        keep_delivery <- is_anchor
      }

      keep_sf <- is_control | (is_delivery & keep_delivery)
      design_sf <- design[keep_sf, , drop = FALSE]

      if (nrow(design_sf) < 4) {
        message("NOTE: phase=anchor size-factor set <4 samples for ", dataset_id,
                "; falling back to ALL samples for size factors.")
        design_sf <- design
      } else {
        message("Size factors computed on phase=anchor for ", dataset_id,
                " with ALL controls + phase-filtered delivery (n=", nrow(design_sf), ").")
      }
    }
  }

  design_sf <- design_sf %>% arrange(match(sample_id, colnames(raw_mat)))
  mat_sf <- raw_mat[, design_sf$sample_id, drop = FALSE]
  coldata_sf <- design_sf %>% column_to_rownames("sample_id")

  dds0 <- DESeqDataSetFromMatrix(countData = mat_sf, colData = coldata_sf, design = ~ 1)
  dds0 <- dds0[rowSums(counts(dds0)) >= 10, ]
  dds0 <- estimateSizeFactors(dds0)

  list(
    sf = sizeFactors(dds0),
    norm_counts = counts(dds0, normalized = TRUE),
    used_samples = colnames(dds0)
  )
}

# -------------------------
# STEP 5B: run phase-specific split contrasts using phase-specific size factors
# -------------------------
run_deseq2_contrasts <- function(dataset_id, raw_mat, sf_vec,
                                 phase = c("anchor","calibration"),
                                 anchor_max = 24,
                                 calib_lower_exc = 24,
                                 calib_upper_exc = 72) {
  phase <- match.arg(phase)

  design_dir <- file.path(split_design_root, paste0(dataset_id, "_design"))
  if (!dir.exists(design_dir)) {
    message("Skip: missing split design folder: ", design_dir)
    return(invisible(list(
      dataset_id = dataset_id, phase = phase,
      n_design_total = 0, n_written = 0, n_skipped = 0, skip_reasons = "missing_split_design_folder"
    )))
  }

  design_files_all <- list.files(design_dir, pattern = "\\.tsv$", full.names = TRUE)
  design_files_all <- design_files_all[sapply(design_files_all, is_real_split_design)]
  if (length(design_files_all) == 0) {
    message("Skip: no per-contrast design TSVs found in: ", design_dir)
    return(invisible(list(
      dataset_id = dataset_id, phase = phase,
      n_design_total = 0, n_written = 0, n_skipped = 0, skip_reasons = "no_split_designs_found"
    )))
  }

  # Filter split designs by __H= in filename
  Hs_all <- vapply(design_files_all, parse_H_from_name, numeric(1))
  phases_all <- vapply(
    Hs_all, phase_from_H, character(1),
    anchor_max = anchor_max,
    calib_lower_exc = calib_lower_exc,
    calib_upper_exc = calib_upper_exc
  )

  keep <- phases_all == phase
  design_files <- design_files_all[keep]

  if (length(design_files) == 0) {
    message("No split designs for phase=", phase, " in ", dataset_id, " (based on __H= in filename).")
    return(invisible(list(
      dataset_id = dataset_id, phase = phase,
      n_design_total = 0, n_written = 0, n_skipped = 0, skip_reasons = "no_designs_in_phase"
    )))
  }

  out_root <- file.path(out_base, dataset_id, "deseq2_contrasts", phase)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  written <- character(0)

  # tracking counters
  n_design_total <- length(design_files)
  skip_reasons <- list()
  skip_inc <- function(reason) {
    skip_reasons[[reason]] <<- (skip_reasons[[reason]] %||% 0) + 1
  }

  for (df_path in design_files) {
    d <- read_table_robust(df_path)

    if (!all(c("sample_id","condition_simple") %in% names(d))) {
      skip_inc("missing_required_columns")
      message("Skip design (missing sample_id/condition_simple): ", basename(df_path))
      next
    }

    d <- d %>%
      mutate(
        sample_id = trimws(as.character(sample_id)),
        condition_simple = factor(as.character(condition_simple), levels = c("CONTROL","DELIVERY"))
      ) %>%
      filter(condition_simple %in% c("CONTROL","DELIVERY"))

    # Must contain BOTH groups
    u_conds <- unique(as.character(d$condition_simple))
    if (length(u_conds) < 2) {
      skip_inc("one_group_only_initial")
      message("Skip design (only one group present: ", paste(u_conds, collapse = ","), "): ", basename(df_path))
      next
    }

    tab2 <- table(d$condition_simple)
    if (any(tab2[c("CONTROL","DELIVERY")] < min_n_per_group)) {
      skip_inc("min_n_per_group_failed")
      message("Skip design (<", min_n_per_group, " per group): ", basename(df_path))
      next
    }

    # Filter to samples present in counts AND in size factors
    samples <- intersect(d$sample_id, intersect(colnames(raw_mat), names(sf_vec)))
    if (length(samples) < 4) {
      skip_inc("too_few_samples_counts_sf")
      message("Skip design (too few samples found in counts/sf): ", basename(df_path))
      next
    }

    d2 <- d %>%
      filter(sample_id %in% samples) %>%
      arrange(match(sample_id, colnames(raw_mat)))

    # After filtering, still need both groups
    u2 <- unique(as.character(d2$condition_simple))
    if (length(u2) < 2) {
      skip_inc("lost_group_after_align")
      message("Skip design (lost one group after aligning to counts/sf): ", basename(df_path))
      next
    }

    counts_sub <- raw_mat[, d2$sample_id, drop = FALSE]
    assert_integer_counts(counts_sub, file_label = paste0(dataset_id, " / ", basename(df_path), " counts_sub"))
    storage.mode(counts_sub) <- "integer"

    coldata <- d2 %>%
      select(sample_id, condition_simple) %>%
      column_to_rownames("sample_id")

    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ condition_simple)
    dds <- dds[rowSums(counts(dds)) >= gene_min_total, ]

    if (nrow(dds) < 50) {
      skip_inc("too_few_genes_after_filter")
      message("Skip design (too few genes after filtering): ", basename(df_path))
      next
    }

    # Ensure size factors exist for all samples
    if (!all(colnames(dds) %in% names(sf_vec))) {
      missing_sf <- setdiff(colnames(dds), names(sf_vec))
      stop("Missing size factors for samples: ", paste(missing_sf, collapse = ", "),
           "\nDataset: ", dataset_id, "\nDesign file: ", df_path)
    }
    if (any(is.na(sf_vec[colnames(dds)]))) {
      stop("NA size factors detected after alignment for ", dataset_id, " / ", basename(df_path))
    }

    sizeFactors(dds) <- sf_vec[colnames(dds)]

    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition_simple", "DELIVERY", "CONTROL"))

    res_df_full <- as.data.frame(res) %>%
      rownames_to_column("gene_id")

    res_df_workflow <- res_df_full %>%
      transmute(
        gene   = gene_id,
        log2FC = log2FoldChange,
        SE     = lfcSE,
        stat   = stat,
        pval   = pvalue,
        FDR    = padj
      )

    stem <- tools::file_path_sans_ext(basename(df_path))
    out_file_full <- file.path(out_root, paste0(dataset_id, "__", stem, "__DE_full.tsv"))
    out_file_wf   <- file.path(out_root, paste0(dataset_id, "__", stem, "__DE_workflow.tsv"))

    write_tsv(res_df_full, out_file_full)
    write_tsv(res_df_workflow, out_file_wf)

    written <- c(written, out_file_wf)
    message("Wrote: ", out_file_wf)
  }

  index_file <- file.path(out_root, paste0(dataset_id, "__contrast_index.tsv"))
  if (length(written) > 0) {
    tibble(result_file = written) %>% write_tsv(index_file)
    message("Index: ", index_file)
  } else {
    message("No contrasts written for ", dataset_id, " phase=", phase)
  }

  skip_str <- if (length(skip_reasons) == 0) {
    ""
  } else {
    paste(paste0(names(skip_reasons), "=", unlist(skip_reasons)), collapse = ";")
  }

  return(invisible(list(
    dataset_id = dataset_id,
    phase = phase,
    n_design_total = n_design_total,
    n_written = length(written),
    n_skipped = n_design_total - length(written),
    skip_reasons = skip_str,
    index_file = index_file
  )))
}

# -------------------------
# PER-DATASET DRIVER
# -------------------------
run_one_dataset <- function(dataset_id) {
  message("\n==============================")
  message("DATASET: ", dataset_id)
  message("==============================")

  # Pre-count split design metadata (for tracker visibility)
  split_counts <- count_split_designs_by_phase(
    dataset_id,
    anchor_max = anchor_hours_max,
    calib_lower_exc = calib_hours_lower_exc,
    calib_upper_exc = calib_hours_upper_exc
  )

  raw_counts_path <- file.path(
    project_root, "03_counts", dataset_id,
    "featurecounts", "validation", "gene_counts_clean.tsv"
  )
  if (!file.exists(raw_counts_path)) {
    message("Skip: missing raw counts: ", raw_counts_path)

    # still track that we skipped dataset
    track_add(list(
      dataset_id = dataset_id,
      phase = "dataset",
      status = "skipped_missing_counts",
      split_total = split_counts$total,
      split_anchor = split_counts$anchor,
      split_calibration = split_counts$calibration,
      split_other = split_counts$other,
      split_missingH = split_counts$missingH
    ))
    return(invisible(NULL))
  }

  full_design_path <- file.path(
    project_root, "00_metadata", "verified_metadata",
    paste0(dataset_id, "_design.tsv")
  )

  # load raw counts
  raw_df  <- read_table_robust(raw_counts_path)
  raw_mat <- make_counts_matrix(raw_df)

  # validate raw integer counts
  assert_integer_counts(raw_mat, file_label = paste0(dataset_id, " raw_mat"))
  storage.mode(raw_mat) <- "integer"

  # ---- ANCHOR phase ----
  if (run_anchor) {
    sf_res_a <- compute_size_factors(
      dataset_id, raw_mat, full_design_path,
      phase = "anchor",
      anchor_max = anchor_hours_max,
      calib_lower_exc = calib_hours_lower_exc,
      calib_upper_exc = calib_hours_upper_exc
    )
    sf_vec_a <- sf_res_a$sf

    norm_out_dir_a <- file.path(project_root, "04_de", dataset_id, "normalized", "anchor")
    dir.create(norm_out_dir_a, showWarnings = FALSE, recursive = TRUE)

    write_tsv(
      tibble(sample_id = names(sf_vec_a), size_factor = as.numeric(sf_vec_a)),
      file.path(norm_out_dir_a, "size_factors.tsv")
    )

    write_tsv(
      as.data.frame(sf_res_a$norm_counts) %>% rownames_to_column("gene_id"),
      file.path(norm_out_dir_a, "gene_counts_normalized.tsv")
    )

    raw_mat_a <- raw_mat[, intersect(colnames(raw_mat), names(sf_vec_a)), drop = FALSE]

    stats_a <- run_deseq2_contrasts(
      dataset_id, raw_mat_a, sf_vec_a,
      phase = "anchor",
      anchor_max = anchor_hours_max,
      calib_lower_exc = calib_hours_lower_exc,
      calib_upper_exc = calib_hours_upper_exc
    )

    track_add(list(
      dataset_id = dataset_id,
      phase = "anchor",
      status = "ran",
      counts_samples_n = ncol(raw_mat),
      sf_used_samples_n = length(sf_res_a$used_samples),
      split_total = split_counts$total,
      split_anchor = split_counts$anchor,
      split_calibration = split_counts$calibration,
      split_other = split_counts$other,
      split_missingH = split_counts$missingH,
      designs_total_in_phase = stats_a$n_design_total,
      contrasts_written = stats_a$n_written,
      contrasts_skipped = stats_a$n_skipped,
      skip_reasons = stats_a$skip_reasons
    ))
  }

  # ---- CALIBRATION phase ----
  if (run_calibration) {
    sf_res_c <- compute_size_factors(
      dataset_id, raw_mat, full_design_path,
      phase = "calibration",
      anchor_max = anchor_hours_max,
      calib_lower_exc = calib_hours_lower_exc,
      calib_upper_exc = calib_hours_upper_exc
    )
    sf_vec_c <- sf_res_c$sf

    norm_out_dir_c <- file.path(project_root, "04_de", dataset_id, "normalized", "calibration")
    dir.create(norm_out_dir_c, showWarnings = FALSE, recursive = TRUE)

    write_tsv(
      tibble(sample_id = names(sf_vec_c), size_factor = as.numeric(sf_vec_c)),
      file.path(norm_out_dir_c, "size_factors.tsv")
    )

    write_tsv(
      as.data.frame(sf_res_c$norm_counts) %>% rownames_to_column("gene_id"),
      file.path(norm_out_dir_c, "gene_counts_normalized.tsv")
    )

    raw_mat_c <- raw_mat[, intersect(colnames(raw_mat), names(sf_vec_c)), drop = FALSE]

    stats_c <- run_deseq2_contrasts(
      dataset_id, raw_mat_c, sf_vec_c,
      phase = "calibration",
      anchor_max = anchor_hours_max,
      calib_lower_exc = calib_hours_lower_exc,
      calib_upper_exc = calib_hours_upper_exc
    )

    track_add(list(
      dataset_id = dataset_id,
      phase = "calibration",
      status = "ran",
      counts_samples_n = ncol(raw_mat),
      sf_used_samples_n = length(sf_res_c$used_samples),
      split_total = split_counts$total,
      split_anchor = split_counts$anchor,
      split_calibration = split_counts$calibration,
      split_other = split_counts$other,
      split_missingH = split_counts$missingH,
      designs_total_in_phase = stats_c$n_design_total,
      contrasts_written = stats_c$n_written,
      contrasts_skipped = stats_c$n_skipped,
      skip_reasons = stats_c$skip_reasons
    ))
  }

  invisible(NULL)
}

# -------------------------
# MAIN
# -------------------------
if (!dir.exists(de_root)) {
  stop("DE root not found: ", de_root,
       "\nExpected: <project_root>/04_de to exist and contain GSE* subfolders.")
}

dataset_ids <- list.dirs(de_root, recursive = FALSE, full.names = FALSE)
dataset_ids <- dataset_ids[grepl("^GSE", dataset_ids)]

if (!is.na(dataset_only)) {
  if (!(dataset_only %in% dataset_ids)) {
    stop("Requested dataset_id not found under 04_de/: ", dataset_only,
         "\nAvailable: ", paste(dataset_ids, collapse = ", "))
  }
  dataset_ids <- dataset_only
}

print(dataset_ids)

for (ds in dataset_ids) run_one_dataset(ds)

track_write(out_base)

message("\nAll Step 5 (phase-specific) normalization + DE contrasts finished.")
message("Normalized outputs under: <project_root>/04_de/<DATASET>/normalized/(anchor|calibration)/")
message("DE outputs under:         ", out_base, "/<DATASET>/deseq2_contrasts/(anchor|calibration)/")
message("Tracker under:            ", out_base, "/_tracker/step5_phase_tracker.tsv")