#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 06A — Core gene sets (LOCKED mouse datasets)
# - Builds core gene set separately for:
#     (1) ANCHOR phase      (H <= 24)    -> Step 5 output: .../anchor/
#     (2) CALIBRATION phase (24 < H <72) -> Step 5 output: .../calibration/
#
# Reads Step 5 DE tables (workflow schema preferred):
#   gene, log2FC, SE, stat, pval, FDR
#
# Also supports DESeq2 "full" tables if present:
#   gene/log2FoldChange/padj (mapped -> log2FC/FDR)
#
# Outputs:
#   <project_root>/05_score/anchors/core_gene_set.tsv
#   <project_root>/05_score/calibration/core_gene_set.tsv
#   plus contrast-count audits per phase
#
# NOTE:
#   This script does NOT read metadata/design files directly.
#   It only reads Step 5 DE outputs under:
#     <project_root>/04_de/comparison/
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

# -------------------------
# USER SETTINGS
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

de_comp_root <- file.path(project_root, "04_de", "comparison")

# Locked mouse datasets
LOCKED_DATASETS_MOUSE <- c(
  "GSE39129",
  "GSE167521",
  "GSE264344",
  "GSE279372",
  "GSE279744",
  "GSE262515"
)

# Thresholds
fdr_cutoff <- 0.05
lfc_cutoff <- 1.0

# Gene must pass in >=50% contrasts within a dataset
within_dataset_support_cutoff <- 0.50

# Gene must be supported by >= ceil(K_present * 2/3) datasets
consensus_k_fraction <- 2 / 3

# Direction control:
#   "both": abs(log2FC) >= lfc_cutoff
#   "up":   log2FC >= lfc_cutoff
#   "down": log2FC <= -lfc_cutoff
direction <- "both"

# Which phases to build
PHASES <- c("anchor", "calibration")

# -------------------------
# Helpers
# -------------------------
read_de_table <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  coln <- names(df)

  gene_col <- dplyr::case_when(
    "gene" %in% coln    ~ "gene",
    "gene_id" %in% coln ~ "gene_id",
    TRUE                ~ NA_character_
  )

  if (is.na(gene_col)) {
    stop(
      "Unrecognized DE table schema in: ", path, "\n",
      "No gene column found. Columns: ", paste(coln, collapse = ", ")
    )
  }

  # Workflow schema
  if (all(c("log2FC", "FDR") %in% coln)) {
    return(
      df %>%
        transmute(
          gene_id = as.character(.data[[gene_col]]),
          log2FC = suppressWarnings(as.numeric(log2FC)),
          FDR = suppressWarnings(as.numeric(FDR))
        )
    )
  }

  # DESeq2 full schema
  if (all(c("log2FoldChange", "padj") %in% coln)) {
    return(
      df %>%
        transmute(
          gene_id = as.character(.data[[gene_col]]),
          log2FC = suppressWarnings(as.numeric(log2FoldChange)),
          FDR = suppressWarnings(as.numeric(padj))
        )
    )
  }

  stop(
    "Unrecognized DE table schema in: ", path, "\n",
    "Columns found: ", paste(coln, collapse = ", "), "\n",
    "Expected workflow schema with log2FC/FDR or full DESeq2 schema with log2FoldChange/padj."
  )
}

get_dataset_from_path <- function(path) {
  # Expected:
  # .../04_de/comparison/<DATASET>/deseq2_contrasts/<phase>/<file>
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/")[[1]]
  i <- which(parts == "comparison")

  if (length(i) == 1 && length(parts) >= i + 1) {
    return(as.character(parts[i + 1]))
  }

  NA_character_
}

pass_by_direction <- function(log2FC, lfc_cutoff, direction) {
  if (direction == "both") {
    return(!is.na(log2FC) & abs(log2FC) >= lfc_cutoff)
  }

  if (direction == "up") {
    return(!is.na(log2FC) & log2FC >= lfc_cutoff)
  }

  if (direction == "down") {
    return(!is.na(log2FC) & log2FC <= -lfc_cutoff)
  }

  stop("direction must be one of: both, up, down")
}

phase_out_dir <- function(phase) {
  if (phase == "anchor") return("anchors")
  if (phase == "calibration") return("calibration")
  stop("Bad phase: ", phase)
}

# -------------------------
# Build core gene set for one phase
# -------------------------
build_core_for_phase <- function(phase) {
  stopifnot(phase %in% c("anchor", "calibration"))

  out_root <- file.path(project_root, "05_score", phase_out_dir(phase))
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  # Collect Step 5 DE files from this phase
  all_files <- list.files(
    de_comp_root,
    pattern = "__(DE_workflow|DE_full)\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  phase_files <- all_files[
    str_detect(
      all_files,
      paste0("deseq2_contrasts[\\\\/]+", phase, "[\\\\/]")
    )
  ]

  if (length(phase_files) == 0) {
    message("PHASE=", phase, ": no DE tables found. Skipping.")
    return(invisible(NULL))
  }

  meta <- tibble(
    file = as.character(phase_files),
    dataset_id = purrr::map_chr(phase_files, get_dataset_from_path)
  ) %>%
    mutate(
      dataset_id = as.character(dataset_id),
      is_human = str_detect(dataset_id, "_HUMAN$")
    ) %>%
    filter(
      !is.na(dataset_id),
      !is_human,
      dataset_id %in% LOCKED_DATASETS_MOUSE
    )

  if (nrow(meta) == 0) {
    message(
      "PHASE=", phase,
      ": found files, but none match locked mouse datasets: ",
      paste(LOCKED_DATASETS_MOUSE, collapse = ", "),
      ". Skipping."
    )
    return(invisible(NULL))
  }

  # -------------------------
  # Audit: contrast counts
  # PATCHED:
  #   force dataset_id to character
  #   force dplyr::count to avoid count() masking/list-column issues
  # -------------------------
  meta <- meta %>%
    mutate(dataset_id = as.character(dataset_id))

  by_ds <- meta %>%
    dplyr::count(dataset_id, name = "n_contrasts") %>%
    dplyr::arrange(dplyr::desc(n_contrasts))

  overall <- tibble(
    phase = phase,
    direction = direction,
    datasets_present = dplyr::n_distinct(as.character(meta$dataset_id)),
    total_contrasts = nrow(meta)
  )

  write_tsv(by_ds, file.path(out_root, "contrast_counts_by_dataset.tsv"))
  write_tsv(overall, file.path(out_root, "contrast_counts_overall.tsv"))

  datasets_present <- sort(unique(as.character(meta$dataset_id)))
  K_present <- length(datasets_present)

  if (K_present < 2) {
    message(
      "PHASE=", phase,
      ": fewer than 2 locked datasets present (K_present=", K_present, "). Skipping this phase.\n",
      "Present: ", paste(datasets_present, collapse = ", ")
    )
    return(invisible(NULL))
  }

  required_support <- ceiling(K_present * consensus_k_fraction)

  message(
    "\nPHASE=", phase,
    " | direction=", direction,
    " | locked datasets present: ", paste(datasets_present, collapse = ", "),
    " | K_present=", K_present,
    " | core requires >= ", required_support, " supporting datasets"
  )

  # -------------------------
  # Load DE tables and add pass flag
  # -------------------------
  de_long <- meta %>%
    mutate(de = purrr::map(file, read_de_table)) %>%
    tidyr::unnest(de) %>%
    mutate(
      pass = !is.na(FDR) &
        FDR <= fdr_cutoff &
        pass_by_direction(log2FC, lfc_cutoff, direction)
    )

  if (nrow(de_long) == 0) {
    stop("PHASE=", phase, ": no DE rows after reading DE tables.")
  }

  # -------------------------
  # Within-dataset reproducibility
  # Collapse evidence within each dataset first
  # -------------------------
  support_by_dataset <- de_long %>%
    mutate(
      dataset_id = as.character(dataset_id),
      gene_id = as.character(gene_id),
      pass = ifelse(is.na(pass), FALSE, pass)
    ) %>%
    group_by(dataset_id, gene_id) %>%
    summarise(
      n_contrasts_dataset = dplyr::n_distinct(file),
      n_pass = sum(pass, na.rm = TRUE),
      support_in_dataset = n_pass / n_contrasts_dataset,
      dataset_support_flag = support_in_dataset >= within_dataset_support_cutoff,
      .groups = "drop"
    )

  # -------------------------
  # Across-dataset consensus
  # -------------------------
  core_support <- support_by_dataset %>%
    group_by(gene_id) %>%
    summarise(
      support_datasets = sum(dataset_support_flag, na.rm = TRUE),
      support_fraction = support_datasets / K_present,
      datasets_supporting = paste(
        sort(unique(dataset_id[dataset_support_flag])),
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(support_datasets), desc(support_fraction), gene_id)

  core_gene_set <- core_support %>%
    filter(support_datasets >= required_support) %>%
    arrange(desc(support_datasets), desc(support_fraction), gene_id)

  if (nrow(core_gene_set) == 0) {
    stop(
      "PHASE=", phase, ": core gene set is EMPTY under thresholds.\n",
      "Try lowering lfc_cutoff or within_dataset_support_cutoff.\n",
      "Current: FDR <= ", fdr_cutoff,
      ", direction = ", direction,
      ", |log2FC| cutoff = ", lfc_cutoff,
      ", within_dataset_support_cutoff = ", within_dataset_support_cutoff,
      ", across >= ceil(K_present * ", consensus_k_fraction, ")"
    )
  }

  out_path <- file.path(out_root, "core_gene_set.tsv")
  write_tsv(core_gene_set, out_path)

  # Extra useful audit outputs
  write_tsv(
    support_by_dataset,
    file.path(out_root, "support_by_dataset.tsv")
  )

  write_tsv(
    core_support,
    file.path(out_root, "core_support_all_genes.tsv")
  )

  message(
    "PHASE=", phase,
    ": wrote core gene set: ", out_path,
    " (n=", nrow(core_gene_set), ")"
  )

  invisible(NULL)
}

# -------------------------
# MAIN
# -------------------------
if (!dir.exists(de_comp_root)) {
  stop(
    "DE comparison root not found: ", de_comp_root,
    "\nExpected: <project_root>/04_de/comparison to exist."
  )
}

for (ph in PHASES) {
  build_core_for_phase(ph)
}

message(
  "\nDONE: Step 06A core gene sets built for phases: ",
  paste(PHASES, collapse = ", ")
)