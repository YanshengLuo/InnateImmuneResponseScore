#!/usr/bin/env Rscript

# Shared helpers for the publication extra analyses. These functions only read
# existing IMRS outputs and write into 05_score/publication_extra_analyses.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
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
AUDIT_RESULTS_DIR <- normalizePath(Sys.getenv("IMRS_LAYER2_AUDIT_RESULTS_DIR",
                                              file.path(PROJECT_ROOT, "audit", "results")),
                                   winslash = "/", mustWork = FALSE)
EXTRA_ROOT <- normalizePath(Sys.getenv("IMRS_LAYER2_EXTRA_ROOT",
                                       file.path(PROJECT_ROOT, "05_score",
                                                 "publication_extra_analyses")),
                            winslash = "/", mustWork = FALSE)
EXTRA_SCRIPT_DIR <- file.path(EXTRA_ROOT, "scripts")
EXTRA_RESULTS_DIR <- file.path(EXTRA_ROOT, "results")
EXTRA_FIGURES_DIR <- file.path(EXTRA_ROOT, "figures")
EXTRA_LOG_DIR <- file.path(EXTRA_ROOT, "logs")

STEP09_EVAL_PATH <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval",
                              "step09_split_eval.tsv")
STEP09_SAMPLE_PATH <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval",
                                "step09_split_sample_level.tsv")
ROLE_TABLE_PATH <- file.path(AUDIT_RESULTS_DIR, "manuscript_dataset_role_table.tsv")

INTERPRETATION_ORDER <- c(
  "strict_anchor",
  "primary_acute_validation",
  "extended_validation",
  "secondary_support_not_primary",
  "excluded_or_unclear"
)

INTERPRETATION_LABELS <- c(
  # Legacy identifier retained for figure compatibility; the role table now
  # assigns all five production anchors to this principal group.
  strict_anchor = "Locked anchor (production)",
  primary_acute_validation = "Primary acute validation",
  extended_validation = "Extended validation",
  secondary_support_not_primary = "Secondary support",
  excluded_or_unclear = "Excluded/unclear"
)

INTERPRETATION_PALETTE <- c(
  strict_anchor = "#2F6B9A",
  primary_acute_validation = "#C84B31",
  extended_validation = "#7B5EA7",
  secondary_support_not_primary = "#59A14F",
  excluded_or_unclear = "#888888"
)

ensure_extra_dirs <- function() {
  invisible(lapply(
    c(EXTRA_ROOT, EXTRA_SCRIPT_DIR, EXTRA_RESULTS_DIR, EXTRA_FIGURES_DIR,
      EXTRA_LOG_DIR),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  ))
}

log_extra <- function(...) {
  ensure_extra_dirs()
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 paste(..., collapse = ""))
  message(line)
  cat(line, "\n", file = file.path(EXTRA_LOG_DIR, "publication_extra_run.log"),
      append = TRUE)
}

read_tsv_if <- function(path) {
  if (!file.exists(path)) return(tibble())
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(x, path, na = "NA")
  invisible(path)
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, path, na = "NA")
  invisible(path)
}

save_plot_pair <- function(plot_obj, stem, width = 8, height = 5, dpi = 300) {
  dir.create(EXTRA_FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(EXTRA_FIGURES_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(EXTRA_FIGURES_DIR, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(pdf_path, plot_obj, width = width, height = height,
                  device = grDevices::cairo_pdf)
  invisible(c(png_path, pdf_path))
}

theme_extra <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
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

safe_quantile <- function(x, probs) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

safe_prop <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

cohens_d_two_group <- function(delivery, control) {
  delivery <- suppressWarnings(as.numeric(delivery))
  control <- suppressWarnings(as.numeric(control))
  delivery <- delivery[is.finite(delivery)]
  control <- control[is.finite(control)]
  if (length(delivery) < 2 || length(control) < 2) return(NA_real_)
  pooled <- sqrt(((length(delivery) - 1) * stats::var(delivery) +
                    (length(control) - 1) * stats::var(control)) /
                   (length(delivery) + length(control) - 2))
  if (!is.finite(pooled) || pooled <= 0) return(NA_real_)
  (mean(delivery) - mean(control)) / pooled
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

strip_ensembl_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

first_present_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

collapse_unique_chr <- function(x, sep = ";") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

interpretation_from_role <- function(manuscript_role) {
  secondary_role <- paste0("cali", "bration")
  dplyr::case_when(
    manuscript_role %in% c("strict_anchor", "anchor") ~ "strict_anchor",
    manuscript_role %in% c("external_acute", "primary_acute_validation") ~
      "primary_acute_validation",
    manuscript_role %in% c("external_extended", "extended_validation") ~
      "extended_validation",
    manuscript_role %in% c(secondary_role, "secondary_support") ~
      "secondary_support_not_primary",
    TRUE ~ "excluded_or_unclear"
  )
}

load_role_table <- function() {
  role_tbl <- read_tsv_if(ROLE_TABLE_PATH)
  if (nrow(role_tbl) == 0) return(tibble())
  role_tbl %>%
    mutate(
      pass = as.logical(pass),
      time_h = suppressWarnings(as.numeric(time_h)),
      manuscript_role = as.character(manuscript_role),
      manuscript_interpretation_group = interpretation_from_role(manuscript_role),
      manuscript_interpretation_label =
        unname(INTERPRETATION_LABELS[manuscript_interpretation_group])
    ) %>%
    distinct(gse_id, split_id, .keep_all = TRUE)
}

load_eval_with_roles <- function(required = TRUE) {
  if (!file.exists(STEP09_EVAL_PATH)) {
    if (required) stop("Missing Step 09 eval table: ", STEP09_EVAL_PATH)
    return(tibble())
  }
  eval_tbl <- readr::read_tsv(STEP09_EVAL_PATH, show_col_types = FALSE,
                              progress = FALSE) %>%
    mutate(
      pass = as.logical(pass),
      time_h = suppressWarnings(as.numeric(time_h)),
      delta_mean_imrs_z = suppressWarnings(as.numeric(delta_mean_imrs_z)),
      auc_imrs_z = suppressWarnings(as.numeric(auc_imrs_z))
    )

  role_tbl <- load_role_table()
  if (nrow(role_tbl) > 0) {
    eval_tbl <- eval_tbl %>%
      left_join(
        role_tbl %>%
          select(gse_id, split_id, manuscript_role, manuscript_claim_group,
                 manuscript_interpretation_group,
                 manuscript_interpretation_label,
                 inclusion_rationale, limitation),
        by = c("gse_id", "split_id")
      )
  }

  if (!("manuscript_interpretation_group" %in% names(eval_tbl))) {
    eval_tbl$manuscript_interpretation_group <- "excluded_or_unclear"
  }
  eval_tbl %>%
    mutate(
      manuscript_interpretation_group = ifelse(
        is.na(manuscript_interpretation_group),
        "excluded_or_unclear",
        manuscript_interpretation_group
      ),
      manuscript_interpretation_group = factor(
        manuscript_interpretation_group,
        levels = INTERPRETATION_ORDER
      ),
      manuscript_interpretation_label =
        unname(INTERPRETATION_LABELS[as.character(manuscript_interpretation_group)])
    )
}

load_sample_with_roles <- function(required = TRUE) {
  if (!file.exists(STEP09_SAMPLE_PATH)) {
    if (required) stop("Missing Step 09 sample-level table: ", STEP09_SAMPLE_PATH)
    return(tibble())
  }
  sample_tbl <- readr::read_tsv(STEP09_SAMPLE_PATH, show_col_types = FALSE,
                                progress = FALSE) %>%
    mutate(
      time_h = suppressWarnings(as.numeric(time_h)),
      imrs_z = suppressWarnings(as.numeric(imrs_z)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    )
  eval_roles <- load_eval_with_roles(required = required) %>%
    select(gse_id, split_id, split_path, manuscript_role,
           manuscript_interpretation_group, manuscript_interpretation_label,
           n_controls, n_delivery) %>%
    distinct(gse_id, split_id, .keep_all = TRUE)
  sample_tbl %>%
    left_join(eval_roles, by = c("gse_id", "split_id")) %>%
    mutate(
      manuscript_interpretation_group = factor(
        ifelse(is.na(manuscript_interpretation_group),
               "excluded_or_unclear",
               as.character(manuscript_interpretation_group)),
        levels = INTERPRETATION_ORDER
      ),
      manuscript_interpretation_label =
        unname(INTERPRETATION_LABELS[as.character(manuscript_interpretation_group)])
    )
}

read_table_robust <- function(path, guess_max = 100000) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  readers <- list(
    function() readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                               guess_max = guess_max),
    function() readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
                               guess_max = guess_max),
    function() readr::read_table(path, show_col_types = FALSE, progress = FALSE,
                                 guess_max = guess_max)
  )
  for (reader in readers) {
    out <- tryCatch(reader(), error = function(e) NULL)
    if (!is.null(out) && ncol(out) >= 1) return(out)
  }
  stop("Could not read table with TSV, CSV, or whitespace parser: ", path)
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
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
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

prepare_dataset_expression <- function(dataset_id, min_controls = 3L) {
  counts_path <- find_counts_file(dataset_id)
  design_path <- find_scoring_design(dataset_id)
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
    filter(sample_id %in% colnames(counts),
           condition_simple %in% c("CONTROL", "DELIVERY")) %>%
    distinct(sample_id, .keep_all = TRUE)
  if (nrow(design) == 0) stop("No overlapping scoring samples for ", dataset_id)
  ctrl_ids <- design$sample_id[design$condition_simple == "CONTROL"]
  if (length(ctrl_ids) < min_controls) {
    stop("Controls < ", min_controls, " for ", dataset_id, " (n=",
         length(ctrl_ids), ")")
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

load_gene_symbol_map <- function() {
  path <- Sys.getenv("IMRS_PORTED_GENE_SYMBOL_MAP",
                     file.path(PROJECT_ROOT, "05_score", "anchors",
                               "gene_symbol_mapping.tsv"))
  if (!file.exists(path)) {
    return(tibble(gene_id = character(), gene_symbol = character()))
  }
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  gene_col <- first_present_col(df, c("ensembl_gene_id", "gene_id", "gene"))
  symbol_col <- first_present_col(df, c("mgi_symbol", "gene_symbol", "symbol",
                                        "external_gene_name"))
  if (is.na(gene_col) || is.na(symbol_col)) {
    return(tibble(gene_id = character(), gene_symbol = character()))
  }
  df %>%
    transmute(
      gene_id = strip_ensembl_version(.data[[gene_col]]),
      gene_symbol = as.character(.data[[symbol_col]])
    ) %>%
    filter(!is.na(gene_id), nzchar(gene_id), !is.na(gene_symbol),
           nzchar(gene_symbol)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

read_split_table <- function(split_path) {
  split <- read_table_robust(split_path)
  if (!("sample_id" %in% names(split))) stop("Split missing sample_id: ", split_path)
  if (!("condition_simple" %in% names(split))) {
    if ("is_control" %in% names(split)) {
      split <- split %>%
        mutate(condition_simple = ifelse(as.logical(is_control), "CONTROL",
                                         "DELIVERY"))
    } else {
      stop("Split missing condition_simple/is_control: ", split_path)
    }
  }
  split %>%
    mutate(
      sample_id = trimws(as.character(sample_id)),
      condition_simple = toupper(trimws(as.character(condition_simple)))
    )
}

inventory_extra_outputs <- function() {
  files <- list.files(EXTRA_ROOT, recursive = TRUE, full.names = TRUE)
  tibble(path = normalizePath(files, winslash = "/", mustWork = FALSE)) %>%
    mutate(
      file = basename(path),
      relative_path = sub(paste0("^", gsub("\\\\", "/", EXTRA_ROOT), "/?"),
                          "", path),
      extension = tools::file_ext(path),
      size_bytes = file.info(path)$size,
      modified_time = as.character(file.info(path)$mtime)
    ) %>%
    arrange(relative_path)
}

ensure_extra_dirs()
