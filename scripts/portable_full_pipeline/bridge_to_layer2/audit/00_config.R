#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("IMRS_PORTED_OUTPUT_ROOT",
                                         Sys.getenv("IMRS_PROJECT_ROOT", ".")),
                               winslash = "/", mustWork = FALSE)
COUNTS_ROOT <- normalizePath(Sys.getenv("IMRS_PORTED_COUNTS_ROOT",
                                        file.path(PROJECT_ROOT, "03_counts")),
                             winslash = "/", mustWork = FALSE)
SCORING_DESIGN_ROOT <- normalizePath(Sys.getenv("IMRS_PORTED_SCORING_DESIGN_ROOT",
                                                file.path(PROJECT_ROOT, "00_metadata",
                                                          "verified_metadata", "scoring")),
                                     winslash = "/", mustWork = FALSE)
AUDIT_ROOT <- normalizePath(Sys.getenv("IMRS_LAYER2_AUDIT_ROOT",
                                       file.path(PROJECT_ROOT, "audit")),
                            winslash = "/", mustWork = FALSE)
AUDIT_SCRIPTS_DIR <- file.path(AUDIT_ROOT, "scripts")
AUDIT_RESULTS_DIR <- file.path(AUDIT_ROOT, "results")
AUDIT_PLOTS_DIR <- file.path(AUDIT_ROOT, "plots")
AUDIT_REPORT_DIR <- file.path(AUDIT_ROOT, "report")
AUDIT_LOG_DIR <- file.path(AUDIT_ROOT, "logs")
AUDIT_DATASETS_DIR <- file.path(AUDIT_ROOT, "datasets")
AUDIT_LOG_FILE <- file.path(AUDIT_LOG_DIR, "audit_run_log.txt")
AUDIT_STATUS_FILE <- file.path(AUDIT_RESULTS_DIR, "audit_section_status.tsv")

# Strict-3 inputs are used only by supporting sensitivity/ablation analyses
# (including threshold sensitivity and strict-3 leave-one-out checks). They do
# not define, replace, or refit the production five-anchor frozen model.
STRICT_ANCHOR_IDS <- c("GSE264344", "GSE39129", "GSE167521")
KNOWN_CALIBRATION_IDS <- c(
  "GSE279372", "GSE279744", "GSE262515", "GSE262515_tissue",
  "GSE262515_cell_line"
)

DEFAULT_MIN_CONTROLS <- 3L
DEFAULT_GENE_SD_FLOOR <- 0.10
DEFAULT_SCORE_SD_FLOOR <- 0.10
DEFAULT_MIN_COVERAGE <- 0.80
DEFAULT_TOP_GENE_N <- 25L

ensure_audit_dirs <- function() {
  dirs <- c(AUDIT_ROOT, AUDIT_SCRIPTS_DIR, AUDIT_RESULTS_DIR, AUDIT_PLOTS_DIR,
            AUDIT_REPORT_DIR, AUDIT_LOG_DIR, AUDIT_DATASETS_DIR)
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
}

initialize_audit_log <- function() {
  ensure_audit_dirs()
  cat(
    paste0("IMRS audit run started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n"),
    file = AUDIT_LOG_FILE
  )
  status <- tibble(
    section = character(),
    status = character(),
    message = character(),
    timestamp = character()
  )
  write_tsv(status, AUDIT_STATUS_FILE)
}

log_msg <- function(level = "INFO", ..., .sep = "") {
  ensure_audit_dirs()
  msg <- paste(..., sep = .sep)
  line <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    toupper(level), " ", msg
  )
  message(line)
  cat(line, "\n", file = AUDIT_LOG_FILE, append = TRUE)
}

record_section_status <- function(section, status, message = NA_character_) {
  ensure_audit_dirs()
  row <- tibble(
    section = section,
    status = status,
    message = ifelse(is.na(message), NA_character_, as.character(message)),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  if (file.exists(AUDIT_STATUS_FILE)) {
    old <- suppressMessages(read_tsv(AUDIT_STATUS_FILE, show_col_types = FALSE))
    old <- old %>%
      mutate(
        section = as.character(section),
        status = as.character(status),
        message = as.character(message),
        timestamp = as.character(timestamp)
      )
    row <- bind_rows(old, row)
  }
  write_tsv(row, AUDIT_STATUS_FILE)
}

read_table_robust <- function(path, guess_max = 100000) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  readers <- list(
    function() read_tsv(path, show_col_types = FALSE, progress = FALSE,
                        guess_max = guess_max),
    function() read_csv(path, show_col_types = FALSE, progress = FALSE,
                        guess_max = guess_max),
    function() read_table(path, show_col_types = FALSE, progress = FALSE,
                          guess_max = guess_max)
  )
  for (reader in readers) {
    out <- tryCatch(reader(), error = function(e) NULL)
    if (!is.null(out) && ncol(out) >= 1) return(out)
  }
  stop("Could not read table with TSV, CSV, or whitespace parser: ", path)
}

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write_tsv(x, path, na = "NA")
  invisible(path)
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write_csv(x, path, na = "NA")
  invisible(path)
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

safe_prop <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

collapse_unique_chr <- function(x, sep = ";") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

first_present_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

strip_ensembl_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

sanitize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_=-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

extract_gse_id_from_path <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/")[[1]]
  hit <- parts[grepl("^GSE[0-9]+(_[A-Za-z0-9]+)?$", parts)]
  if (length(hit) == 0) return(NA_character_)
  hit[[length(hit)]]
}

extract_gse_id_from_split_id <- function(split_id) {
  out <- str_match(as.character(split_id), "(GSE[0-9]+(?:_[A-Za-z0-9]+)?)")[, 2]
  ifelse(is.na(out), NA_character_, out)
}

extract_H_from_filename <- function(path) {
  b <- basename(path)
  m <- str_match(b, "__H=([^_]+)(?:__|$)")
  val <- suppressWarnings(as.numeric(m[, 2]))
  ifelse(length(val) == 0 || is.na(val), NA_real_, val)
}

assign_audit_role <- function(gse_id, time_h = NA_real_) {
  gse_id <- as.character(gse_id)
  time_h <- suppressWarnings(as.numeric(time_h))
  out <- rep("external", length(gse_id))
  out[gse_id %in% KNOWN_CALIBRATION_IDS] <- "calibration"
  out[gse_id %in% setdiff(STRICT_ANCHOR_IDS, "GSE264344")] <- "anchor"
  is_gse264344 <- gse_id == "GSE264344"
  out[is_gse264344 & (!is.finite(time_h) | time_h <= 24)] <- "anchor"
  out[is_gse264344 & is.finite(time_h) & time_h > 24] <- "external"
  out
}

assign_time_group <- function(time_h) {
  time_h <- suppressWarnings(as.numeric(time_h))
  case_when(
    !is.finite(time_h) ~ "unknown_time",
    time_h <= 24 ~ "acute_main",
    time_h > 24 ~ "extended_transfer",
    TRUE ~ "unknown_time"
  )
}

load_step09_eval <- function(required = TRUE) {
  path <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval", "step09_split_eval.tsv")
  if (!file.exists(path)) {
    if (required) stop("Missing Step 09 eval table: ", path)
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      time_h = suppressWarnings(as.numeric(time_h)),
      audit_role = assign_audit_role(gse_id, time_h),
      time_group = assign_time_group(time_h)
    )
}

load_step09_sample_level <- function(required = TRUE) {
  path <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval",
                    "step09_split_sample_level.tsv")
  if (!file.exists(path)) {
    if (required) stop("Missing Step 09 sample-level table: ", path)
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      time_h = suppressWarnings(as.numeric(time_h)),
      audit_role = assign_audit_role(gse_id, time_h),
      time_group = assign_time_group(time_h)
    )
}

load_step09_summary <- function(required = TRUE) {
  path <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval",
                    "step09_split_summary.tsv")
  if (!file.exists(path)) {
    if (required) stop("Missing Step 09 summary table: ", path)
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

load_gene_symbol_map <- function() {
  path <- Sys.getenv("IMRS_PORTED_GENE_SYMBOL_MAP",
                     file.path(PROJECT_ROOT, "05_score", "anchors",
                               "gene_symbol_mapping.tsv"))
  if (!file.exists(path)) return(tibble(gene_id = character(), gene_symbol = character()))
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  gene_col <- first_present_col(df, c("ensembl_gene_id", "gene_id", "gene"))
  symbol_col <- first_present_col(df, c("mgi_symbol", "gene_symbol", "symbol", "external_gene_name"))
  if (is.na(gene_col) || is.na(symbol_col)) {
    return(tibble(gene_id = character(), gene_symbol = character()))
  }
  df %>%
    transmute(
      gene_id = strip_ensembl_version(.data[[gene_col]]),
      gene_symbol = as.character(.data[[symbol_col]])
    ) %>%
    filter(!is.na(gene_id), nzchar(gene_id), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

load_frozen_weights <- function(path = file.path(PROJECT_ROOT, "05_score", "anchors",
                                                 "gene_weights.tsv")) {
  if (!file.exists(path)) stop("Missing frozen weights: ", path)
  w <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  gene_col <- first_present_col(w, c("gene_id", "gene", "ensembl_gene_id"))
  beta_col <- first_present_col(w, c("beta_meta", "weight", "beta", "log2FC"))
  if (is.na(gene_col)) stop("Weights table has no gene column: ", path)
  if (is.na(beta_col)) stop("Weights table has no beta/weight column: ", path)
  w <- w %>%
    transmute(
      gene_id_versioned = as.character(.data[[gene_col]]),
      gene_id = strip_ensembl_version(.data[[gene_col]]),
      beta_meta = suppressWarnings(as.numeric(.data[[beta_col]])),
      source_weight = if ("weight" %in% names(.)) suppressWarnings(as.numeric(weight)) else beta_meta
    ) %>%
    filter(!is.na(gene_id), nzchar(gene_id), is.finite(beta_meta)) %>%
    group_by(gene_id) %>%
    slice_max(abs(beta_meta), n = 1, with_ties = FALSE) %>%
    ungroup()
  symbols <- load_gene_symbol_map()
  if (nrow(symbols) > 0) w <- left_join(w, symbols, by = "gene_id")
  if (!("gene_symbol" %in% names(w))) w$gene_symbol <- NA_character_
  w %>% arrange(desc(abs(beta_meta)))
}

find_counts_file <- function(dataset_id) {
  # PORTING NOTE: Both manuscript GSE262515 arm identifiers use the same prepared count matrix.
  count_dataset_id <- if (dataset_id %in% c("GSE262515_tissue", "GSE262515_cell_line")) {
    "GSE262515"
  } else {
    dataset_id
  }
  root <- file.path(COUNTS_ROOT, count_dataset_id)
  if (!dir.exists(root)) return(NA_character_)
  files <- list.files(root, pattern = "^gene_counts_clean\\.tsv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files <- files[order(!grepl("validation", files, ignore.case = TRUE), files)]
  normalizePath(files[[1]], winslash = "/", mustWork = FALSE)
}

find_scoring_design <- function(dataset_id) {
  candidates <- c(
    file.path(SCORING_DESIGN_ROOT, dataset_id, paste0(dataset_id, "_design.tsv")),
    file.path(SCORING_DESIGN_ROOT, "..",
              paste0(dataset_id, "_design.tsv"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = FALSE)
}

read_counts_matrix <- function(path, strip_versions = TRUE) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (ncol(df) < 2) stop("Counts file has fewer than two columns: ", path)
  gene_id <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- if (strip_versions) strip_ensembl_version(gene_id) else gene_id
  colnames(mat) <- trimws(colnames(mat))
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  mat
}

is_integerish_matrix <- function(mat) {
  all(is.finite(mat)) && all(abs(mat - round(mat)) < 1e-6)
}

estimate_size_factors_imrs <- function(counts_int) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for Step 08-compatible size-factor normalization.")
  }
  DESeq2::estimateSizeFactorsForMatrix(counts_int)
}

score_from_z <- function(z, wvec, design, coverage,
                         score_sd_floor = DEFAULT_SCORE_SD_FLOOR) {
  wvec <- wvec[rownames(z)]
  raw <- as.numeric(crossprod(wvec, z))
  names(raw) <- colnames(z)
  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  ctrl_scores <- raw[ctrl_ids]
  mu_score <- mean(ctrl_scores, na.rm = TRUE)
  sd_score <- sd(ctrl_scores, na.rm = TRUE)
  sd_score_floor <- max(sd_score, score_sd_floor, na.rm = TRUE)
  if (!is.finite(sd_score_floor) || sd_score_floor <= 0) sd_score_floor <- score_sd_floor
  imrs_z <- (raw - mu_score) / sd_score_floor
  tibble(
    sample_id = design$sample_id,
    condition_simple = design$condition_simple,
    imrs_raw = as.numeric(raw[design$sample_id]),
    imrs_z = as.numeric(imrs_z[design$sample_id]),
    coverage = coverage
  )
}

score_dataset_components <- function(dataset_id, weights_df,
                                     counts_path = find_counts_file(dataset_id),
                                     design_path = find_scoring_design(dataset_id),
                                     min_controls = DEFAULT_MIN_CONTROLS,
                                     min_coverage = DEFAULT_MIN_COVERAGE,
                                     gene_sd_floor = DEFAULT_GENE_SD_FLOOR,
                                     score_sd_floor = DEFAULT_SCORE_SD_FLOOR,
                                     enforce_min_coverage = TRUE) {
  if (is.na(counts_path) || !file.exists(counts_path)) {
    stop("Missing counts for dataset ", dataset_id)
  }
  if (is.na(design_path) || !file.exists(design_path)) {
    stop("Missing scoring design for dataset ", dataset_id)
  }
  counts <- read_counts_matrix(counts_path, strip_versions = TRUE)
  design <- read_table_robust(design_path) %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    ) %>%
    filter(sample_id %in% colnames(counts)) %>%
    distinct(sample_id, .keep_all = TRUE)
  if (!all(c("sample_id", "condition_simple") %in% names(design))) {
    stop("Design must contain sample_id and condition_simple: ", design_path)
  }
  if (nrow(design) == 0) stop("No overlapping samples for ", dataset_id)
  design <- design %>% filter(condition_simple %in% c("CONTROL", "DELIVERY"))
  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  if (length(ctrl_ids) < min_controls) {
    stop("Controls < ", min_controls, " for ", dataset_id, " (n=", length(ctrl_ids), ")")
  }
  counts <- counts[, design$sample_id, drop = FALSE]
  if (!is_integerish_matrix(counts)) stop("Counts not integerish for ", dataset_id)
  counts_int <- round(counts)
  storage.mode(counts_int) <- "integer"
  sf <- estimate_size_factors_imrs(counts_int)
  norm <- sweep(counts_int, 2, sf, "/")
  x <- log2(norm + 1)
  w_genes <- unique(weights_df$gene_id)
  overlap <- intersect(rownames(x), w_genes)
  coverage <- length(overlap) / length(w_genes)
  if (enforce_min_coverage && (!is.finite(coverage) || coverage < min_coverage)) {
    stop("Coverage < ", min_coverage, " for ", dataset_id, ": ",
         length(overlap), "/", length(w_genes))
  }
  wsub <- weights_df %>%
    filter(gene_id %in% overlap) %>%
    arrange(match(gene_id, overlap))
  xsub <- x[wsub$gene_id, , drop = FALSE]
  ctrl_mat <- xsub[, ctrl_ids, drop = FALSE]
  mu <- rowMeans(ctrl_mat, na.rm = TRUE)
  sdv <- apply(ctrl_mat, 1, sd, na.rm = TRUE)
  sdv_floor <- pmax(sdv, gene_sd_floor)
  z <- sweep(xsub, 1, mu, "-")
  z <- sweep(z, 1, sdv_floor, "/")
  wvec <- wsub$beta_meta
  names(wvec) <- wsub$gene_id
  score_tbl <- score_from_z(z, wvec, design, coverage, score_sd_floor)
  list(
    dataset_id = dataset_id,
    counts_path = counts_path,
    design_path = design_path,
    design = design,
    z = z,
    wvec = wvec,
    score_tbl = score_tbl,
    coverage = coverage,
    n_overlap = length(overlap),
    n_weights = length(w_genes),
    n_genes_sd_floor = sum(sdv < gene_sd_floor, na.rm = TRUE)
  )
}

prepare_dataset_expression <- function(dataset_id,
                                       counts_path = find_counts_file(dataset_id),
                                       design_path = find_scoring_design(dataset_id),
                                       min_controls = DEFAULT_MIN_CONTROLS) {
  if (is.na(counts_path) || !file.exists(counts_path)) {
    stop("Missing counts for dataset ", dataset_id)
  }
  if (is.na(design_path) || !file.exists(design_path)) {
    stop("Missing scoring design for dataset ", dataset_id)
  }
  counts <- read_counts_matrix(counts_path, strip_versions = TRUE)
  design <- read_table_robust(design_path) %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    ) %>%
    filter(sample_id %in% colnames(counts), condition_simple %in% c("CONTROL", "DELIVERY")) %>%
    distinct(sample_id, .keep_all = TRUE)
  if (nrow(design) == 0) stop("No overlapping scoring samples for ", dataset_id)
  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  if (length(ctrl_ids) < min_controls) {
    stop("Controls < ", min_controls, " for ", dataset_id, " (n=", length(ctrl_ids), ")")
  }
  counts <- counts[, design$sample_id, drop = FALSE]
  if (!is_integerish_matrix(counts)) stop("Counts not integerish for ", dataset_id)
  counts_int <- round(counts)
  storage.mode(counts_int) <- "integer"
  sf <- estimate_size_factors_imrs(counts_int)
  norm <- sweep(counts_int, 2, sf, "/")
  x <- log2(norm + 1)
  list(
    dataset_id = dataset_id,
    counts_path = counts_path,
    design_path = design_path,
    design = design,
    x = x
  )
}

score_prepared_dataset <- function(prep, weights_df,
                                   min_coverage = DEFAULT_MIN_COVERAGE,
                                   gene_sd_floor = DEFAULT_GENE_SD_FLOOR,
                                   score_sd_floor = DEFAULT_SCORE_SD_FLOOR,
                                   enforce_min_coverage = TRUE) {
  x <- prep$x
  design <- prep$design
  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  w_genes <- unique(weights_df$gene_id)
  overlap <- intersect(rownames(x), w_genes)
  coverage <- length(overlap) / length(w_genes)
  if (enforce_min_coverage && (!is.finite(coverage) || coverage < min_coverage)) {
    stop("Coverage < ", min_coverage, " for ", prep$dataset_id, ": ",
         length(overlap), "/", length(w_genes))
  }
  wsub <- weights_df %>%
    filter(gene_id %in% overlap) %>%
    arrange(match(gene_id, overlap))
  xsub <- x[wsub$gene_id, , drop = FALSE]
  ctrl_mat <- xsub[, ctrl_ids, drop = FALSE]
  mu <- rowMeans(ctrl_mat, na.rm = TRUE)
  sdv <- apply(ctrl_mat, 1, sd, na.rm = TRUE)
  sdv_floor <- pmax(sdv, gene_sd_floor)
  z <- sweep(xsub, 1, mu, "-")
  z <- sweep(z, 1, sdv_floor, "/")
  wvec <- wsub$beta_meta
  names(wvec) <- wsub$gene_id
  score_tbl <- score_from_z(z, wvec, design, coverage, score_sd_floor)
  list(
    dataset_id = prep$dataset_id,
    counts_path = prep$counts_path,
    design_path = prep$design_path,
    design = design,
    z = z,
    wvec = wvec,
    score_tbl = score_tbl,
    coverage = coverage,
    n_overlap = length(overlap),
    n_weights = length(w_genes),
    n_genes_sd_floor = sum(sdv < gene_sd_floor, na.rm = TRUE)
  )
}

compute_auc_safe <- function(scores, labels) {
  scores <- suppressWarnings(as.numeric(scores))
  labels <- suppressWarnings(as.integer(labels))
  keep <- is.finite(scores) & !is.na(labels)
  if (sum(keep) == 0 || length(unique(labels[keep])) < 2) return(NA_real_)
  if (!requireNamespace("pROC", quietly = TRUE)) return(NA_real_)
  roc_obj <- tryCatch(
    pROC::roc(response = labels[keep], predictor = scores[keep],
              levels = c(0, 1), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
}

read_split_table <- function(split_path) {
  df <- read_table_robust(split_path)
  if (!("sample_id" %in% names(df))) stop("Split missing sample_id: ", split_path)
  if (!("condition_simple" %in% names(df))) {
    if ("is_control" %in% names(df)) {
      df <- df %>%
        mutate(condition_simple = ifelse(as.logical(is_control), "CONTROL", "DELIVERY"))
    } else {
      stop("Split missing condition_simple/is_control: ", split_path)
    }
  }
  df %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    )
}

compute_eval_time_h <- function(split_path, merged_df) {
  from_name <- extract_H_from_filename(split_path)
  if (is.finite(from_name)) return(from_name)
  if ("time_h" %in% names(merged_df)) {
    del <- suppressWarnings(as.numeric(merged_df$time_h[merged_df$condition_simple == "DELIVERY"]))
    del <- del[is.finite(del)]
    if (length(del) > 0) return(mean(unique(del)))
  }
  NA_real_
}

evaluate_scores_for_split <- function(score_tbl, split_path, gse_id = NA_character_,
                                      score_col = "imrs_z") {
  split <- read_split_table(split_path)
  if (!(score_col %in% names(score_tbl))) stop("Score column missing: ", score_col)
  merged <- split %>%
    select(any_of(c("sample_id", "condition_simple", "tissue", "time_h",
                    "control_label", "contrast_label", "group", "group_raw",
                    "tissue_raw"))) %>%
    inner_join(score_tbl %>% select(sample_id, imrs_raw, imrs_z, coverage),
               by = "sample_id")
  if (nrow(merged) == 0) {
    return(tibble(
      gse_id = gse_id,
      split_id = sub("\\.tsv$", "", basename(split_path)),
      split_path = split_path,
      pass = FALSE,
      fail_reason = "no_scored_samples_overlap"
    ))
  }
  merged <- merged %>%
    mutate(label = ifelse(condition_simple == "DELIVERY", 1L, 0L))
  ctrl <- merged[[score_col]][merged$condition_simple == "CONTROL"]
  del <- merged[[score_col]][merged$condition_simple == "DELIVERY"]
  time_h <- compute_eval_time_h(split_path, merged)
  tibble(
    gse_id = gse_id,
    split_id = sub("\\.tsv$", "", basename(split_path)),
    split_path = split_path,
    pass = length(ctrl) > 0 && length(del) > 0,
    fail_reason = ifelse(length(ctrl) > 0 && length(del) > 0, NA_character_,
                         "missing_control_or_delivery"),
    contrast_label = collapse_unique_chr(merged$contrast_label),
    control_label = collapse_unique_chr(merged$control_label),
    tissue = collapse_unique_chr(merged$tissue),
    time_h = time_h,
    n_split_samples = nrow(split),
    n_scored_samples = nrow(merged),
    n_controls = sum(merged$condition_simple == "CONTROL"),
    n_delivery = sum(merged$condition_simple == "DELIVERY"),
    mean_coverage = safe_mean(merged$coverage),
    control_mean_imrs_z = safe_mean(ctrl),
    delivery_mean_imrs_z = safe_mean(del),
    delta_mean_imrs_z = safe_mean(del) - safe_mean(ctrl),
    control_median_imrs_z = safe_median(ctrl),
    delivery_median_imrs_z = safe_median(del),
    delta_median_imrs_z = safe_median(del) - safe_median(ctrl),
    auc_imrs_z = compute_auc_safe(merged[[score_col]], merged$label)
  )
}

evaluate_scores_for_step09_rows <- function(score_tbl, eval_rows) {
  if (nrow(eval_rows) == 0) return(tibble())
  map_dfr(seq_len(nrow(eval_rows)), function(i) {
    row <- eval_rows[i, ]
    tryCatch(
      evaluate_scores_for_split(score_tbl, row$split_path, row$gse_id),
      error = function(e) {
        log_msg("WARN", "Split evaluation failed for ", row$split_path, ": ",
                conditionMessage(e))
        tibble(
          gse_id = row$gse_id,
          split_id = row$split_id,
          split_path = row$split_path,
          pass = FALSE,
          fail_reason = conditionMessage(e)
        )
      }
    )
  })
}

find_de_tables <- function(dataset_ids = NULL, phases = NULL, prefer_workflow = TRUE) {
  root <- file.path(PROJECT_ROOT, "04_de", "comparison")
  if (!dir.exists(root)) return(tibble())
  files <- list.files(root, pattern = "__(DE_workflow|DE_full)\\.tsv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(tibble())
  meta <- tibble(path = normalizePath(files, winslash = "/", mustWork = FALSE)) %>%
    mutate(
      dataset_id = map_chr(path, extract_gse_id_from_path),
      phase = case_when(
        str_detect(path, "/anchor/") ~ "anchor",
        str_detect(path, "/calibration/") ~ "calibration",
        TRUE ~ NA_character_
      ),
      table_type = ifelse(str_detect(path, "DE_workflow"), "workflow", "full"),
      contrast_id = str_remove(basename(path), "__(DE_workflow|DE_full)\\.tsv$")
    )
  if (!is.null(dataset_ids)) meta <- meta %>% filter(dataset_id %in% dataset_ids)
  if (!is.null(phases)) meta <- meta %>% filter(phase %in% phases)
  if (prefer_workflow) {
    meta <- meta %>%
      group_by(dataset_id, phase, contrast_id) %>%
      arrange(match(table_type, c("workflow", "full")), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()
  }
  meta
}

read_de_table_flexible <- function(path) {
  df <- read_table_robust(path)
  gene_col <- first_present_col(df, c("gene_id", "gene", "gene_symbol",
                                      "symbol", "ensembl_gene_id"))
  lfc_col <- first_present_col(df, c("log2FoldChange", "log2fc", "log2FC",
                                     "logFC", "lfc", "LFC"))
  se_col <- first_present_col(df, c("lfcSE", "SE", "se", "stdErr", "stderr"))
  fdr_col <- first_present_col(df, c("padj", "FDR", "fdr", "qvalue", "q_value",
                                     "adj.P.Val"))
  if (is.na(gene_col) || is.na(lfc_col)) {
    stop("DE table requires gene and log2FC columns: ", path)
  }
  out <- tibble(
    gene_id = strip_ensembl_version(df[[gene_col]]),
    log2FC = suppressWarnings(as.numeric(df[[lfc_col]])),
    SE = if (!is.na(se_col)) suppressWarnings(as.numeric(df[[se_col]])) else NA_real_,
    FDR = if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]])) else NA_real_
  ) %>%
    filter(!is.na(gene_id), nzchar(gene_id), is.finite(log2FC))
  out
}

build_weights_from_de <- function(training_anchor_ids,
                                  min_abs_log2fc = 0.50,
                                  fdr_support = 0.05,
                                  min_up_anchor_count = length(training_anchor_ids),
                                  within_dataset_support = 0.50,
                                  require_up = TRUE) {
  meta <- find_de_tables(dataset_ids = training_anchor_ids, phases = "anchor",
                         prefer_workflow = TRUE)
  if (nrow(meta) == 0) {
    stop("No anchor DE tables found for: ", paste(training_anchor_ids, collapse = ", "))
  }
  de_long <- meta %>%
    mutate(de = map(path, read_de_table_flexible)) %>%
    unnest(de) %>%
    mutate(
      pass_effect = if (require_up) {
        !is.na(log2FC) & log2FC >= min_abs_log2fc
      } else {
        !is.na(log2FC) & abs(log2FC) >= min_abs_log2fc
      },
      pass = pass_effect & (is.na(FDR) | FDR <= fdr_support)
    )
  support <- de_long %>%
    group_by(dataset_id, gene_id) %>%
    summarise(
      n_contrasts_dataset = n_distinct(path),
      n_pass = sum(pass, na.rm = TRUE),
      support_fraction = n_pass / n_contrasts_dataset,
      dataset_support_flag = support_fraction >= within_dataset_support,
      .groups = "drop"
    )
  selected <- support %>%
    group_by(gene_id) %>%
    summarise(
      support_datasets = sum(dataset_support_flag, na.rm = TRUE),
      datasets_supporting = paste(sort(dataset_id[dataset_support_flag]), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(support_datasets >= min_up_anchor_count)
  if (nrow(selected) == 0) {
    stop("No genes selected from training anchors under thresholds.")
  }
  beta <- de_long %>%
    semi_join(selected, by = "gene_id") %>%
    group_by(gene_id) %>%
    summarise(
      beta_meta = {
        ok <- is.finite(SE) & SE > 0
        if (sum(ok) > 0) {
          ww <- 1 / (SE[ok]^2)
          sum(log2FC[ok] * ww) / sum(ww)
        } else {
          mean(log2FC, na.rm = TRUE)
        }
      },
      se_meta = {
        ok <- is.finite(SE) & SE > 0
        if (sum(ok) > 0) sqrt(1 / sum(1 / (SE[ok]^2))) else NA_real_
      },
      I2 = {
        ok <- is.finite(SE) & SE > 0
        if (sum(ok) > 1) {
          ww <- 1 / (SE[ok]^2)
          mu <- sum(log2FC[ok] * ww) / sum(ww)
          q <- sum(ww * (log2FC[ok] - mu)^2)
          max(0, (q - (sum(ok) - 1)) / q)
        } else {
          0
        }
      },
      n_de_rows = n(),
      .groups = "drop"
    ) %>%
    left_join(selected, by = "gene_id") %>%
    filter(is.finite(beta_meta))
  cap_c <- as.numeric(quantile(abs(beta$beta_meta), 0.95, na.rm = TRUE, names = FALSE))
  if (!is.finite(cap_c) || cap_c <= 0) cap_c <- max(abs(beta$beta_meta), na.rm = TRUE)
  beta %>%
    mutate(
      beta_capped = pmax(pmin(beta_meta, cap_c), -cap_c),
      penalty_het = ifelse(is.na(I2), 1, pmax(0, 1 - I2)),
      penalty_final = pmax(penalty_het, 0.2),
      source_weight = beta_capped * penalty_final,
      cap_c = cap_c
    ) %>%
    arrange(desc(abs(beta_meta)))
}

save_plot <- function(plot_obj, path, width = 8, height = 5, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
  invisible(path)
}

theme_audit <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

ensure_audit_dirs()
