#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 06B — Heterogeneity screening (ANCHOR ONLY)
# - Uses Step 5 workflow DE tables: gene, log2FC, SE, stat, pval, FDR
# - Phase: anchor only (H <= 24)
# - Collapses within each dataset first (inverse-variance FE)
# - Then meta across datasets (inverse-variance FE)
# - Computes Q and I2 across DATASETS (not contrasts)
# - Computes for ALL genes, but flags heterogeneity only for core genes
# - Flags only when k_datasets_used >= 3 (stability)
#
# Outputs:
#   <project_root>/05_score/anchors/gene_heterogeneity.tsv
#   plus optional audit: contrast_counts_by_dataset.tsv
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

# Locked mouse datasets (same locked anchor set as Step 06A)
LOCKED_DATASETS_MOUSE <- c("GSE39129", "GSE167521", "GSE264344", "GSE279372", "GSE279744")

# Phase fixed to ANCHOR ONLY
PHASE <- "anchor"

# Heterogeneity flag threshold
I2_flag_threshold <- 0.70

# Numerical safety
se_floor <- 1e-6
max_abs_beta <- 50

# -------------------------
# Helpers
# -------------------------
phase_out_dir <- function(phase) {
  if (phase == "anchor") return("anchors")
  stop("bad phase: ", phase)
}

read_de_workflow_full <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  required <- c("gene", "log2FC", "SE", "stat", "pval", "FDR")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in: ", path, "\nMissing: ", paste(missing, collapse = ", "))
  }

  df %>%
    transmute(
      gene_id = as.character(gene),
      beta    = suppressWarnings(as.numeric(log2FC)),
      se      = suppressWarnings(as.numeric(SE))
    ) %>%
    mutate(
      beta = pmax(pmin(beta,  max_abs_beta), -max_abs_beta),
      se   = ifelse(is.na(se), NA_real_, pmax(se, se_floor))
    )
}

get_dataset_from_path <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/"), "/")[[1]]
  i <- which(parts == "comparison")
  if (length(i) == 1 && length(parts) >= i + 1) return(parts[i + 1])
  NA_character_
}

fe_meta <- function(beta, se) {
  ok <- !is.na(beta) & !is.na(se) & se > 0
  if (!any(ok)) return(list(beta = NA_real_, se = NA_real_, sum_w = 0))

  w <- 1 / (se[ok]^2)
  sum_w <- sum(w)
  beta_m <- sum(w * beta[ok]) / sum_w
  se_m <- 1 / sqrt(sum_w)
  list(beta = beta_m, se = se_m, sum_w = sum_w)
}

# -------------------------
# MAIN
# -------------------------
if (!dir.exists(de_comp_root)) {
  stop("DE comparison root not found: ", de_comp_root,
       "\nExpected: <project_root>/04_de/comparison to exist.")
}

out_root <- file.path(project_root, "05_score", phase_out_dir(PHASE))
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Load core gene set for anchor (for flagging)
core_path <- file.path(out_root, "core_gene_set.tsv")
core_genes <- character(0)
if (file.exists(core_path)) {
  core_df <- read_tsv(core_path, show_col_types = FALSE)
  if ("gene_id" %in% names(core_df)) core_genes <- as.character(core_df$gene_id)
} else {
  message("NOTE: core_gene_set.tsv not found at: ", core_path, "\nHeterogeneity will still be computed for all genes, but no core flagging.")
}

# Collect Step 5 workflow files for ANCHOR
all_files <- list.files(
  de_comp_root,
  pattern = "__DE_workflow\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

phase_files <- all_files[str_detect(all_files, paste0("deseq2_contrasts[\\\\/]+", PHASE, "[\\\\/]"))]

if (length(phase_files) == 0) {
  stop("[6B] PHASE=", PHASE, ": no DE_workflow files found. Did Step 5 write anchor outputs?")
}

meta <- tibble(
  file = phase_files,
  dataset_id = map_chr(phase_files, get_dataset_from_path)
) %>%
  mutate(is_human = str_detect(dataset_id, "_HUMAN$")) %>%
  filter(!is_human, dataset_id %in% LOCKED_DATASETS_MOUSE)

if (nrow(meta) == 0) {
  stop("[6B] PHASE=", PHASE, ": no files match locked mouse datasets.\nLocked: ",
       paste(LOCKED_DATASETS_MOUSE, collapse = ", "))
}

# Audit contrast counts (helpful for debugging dominance)
by_ds <- meta %>%
  count(dataset_id, name = "n_contrasts") %>%
  arrange(desc(n_contrasts), dataset_id)
write_tsv(by_ds, file.path(out_root, "contrast_counts_by_dataset_6B.tsv"))

datasets_present <- sort(unique(meta$dataset_id))
K_present <- length(datasets_present)
if (K_present < 3) {
  stop("[6B] PHASE=", PHASE, ": need >=3 datasets for stable heterogeneity screening.\nPresent: ",
       paste(datasets_present, collapse = ", "))
}

message("\n[6B] PHASE=", PHASE,
        " | locked datasets present: ", paste(datasets_present, collapse = ", "),
        " | K_present=", K_present)

# Long table: dataset_id, contrast_file, gene_id, beta, se
de_long <- meta %>%
  mutate(de = map(file, read_de_workflow_full)) %>%
  tidyr::unnest(de) %>%
  mutate(contrast_file = file)

# Collapse WITHIN dataset first (FE across contrasts)
ds_collapsed <- de_long %>%
  group_by(dataset_id, gene_id) %>%
  summarise(
    n_contrasts = n_distinct(contrast_file),
    tmp = list(fe_meta(beta, se)),
    beta_ds = tmp[[1]]$beta,
    se_ds   = tmp[[1]]$se,
    sum_w_ds = tmp[[1]]$sum_w,
    .groups = "drop"
  ) %>%
  select(-tmp) %>%
  mutate(w_ds = ifelse(!is.na(se_ds) & se_ds > 0, 1/(se_ds^2), NA_real_))

# Meta ACROSS datasets (FE) per gene
gene_meta <- ds_collapsed %>%
  group_by(gene_id) %>%
  summarise(
    datasets_used = paste(sort(unique(dataset_id[!is.na(beta_ds) & !is.na(se_ds) & se_ds > 0])), collapse = ";"),
    k_datasets_used = sum(!is.na(beta_ds) & !is.na(se_ds) & se_ds > 0),
    tmp = list(fe_meta(beta_ds, se_ds)),
    beta_meta = tmp[[1]]$beta,
    se_meta   = tmp[[1]]$se,
    sum_w     = tmp[[1]]$sum_w,
    .groups = "drop"
  ) %>%
  select(-tmp)

# Q and I2 across datasets
gene_qi2 <- ds_collapsed %>%
  inner_join(gene_meta %>% select(gene_id, beta_meta), by = "gene_id") %>%
  filter(!is.na(beta_ds) & !is.na(se_ds) & se_ds > 0) %>%
  group_by(gene_id) %>%
  summarise(
    k = n(),
    Q = sum((1/(se_ds^2)) * (beta_ds - beta_meta)^2),
    df = pmax(k - 1, 0),
    I2 = case_when(
      k < 2 ~ NA_real_,
      Q <= 0 ~ 0,
      TRUE ~ pmax(0, (Q - (k - 1)) / Q)
    ),
    .groups = "drop"
  )

out <- gene_meta %>%
  left_join(gene_qi2, by = "gene_id") %>%
  mutate(
    phase = PHASE,
    is_core = gene_id %in% core_genes,
    heterogeneity_flag = is_core & (k_datasets_used >= 3) & !is.na(I2) & (I2 >= I2_flag_threshold)
  ) %>%
  arrange(desc(is_core), desc(heterogeneity_flag), desc(I2), gene_id)

out_path <- file.path(out_root, "gene_heterogeneity.tsv")
write_tsv(out, out_path)

message("[6B] Wrote: ", out_path,
        " | genes=", nrow(out),
        " | core_genes_loaded=", length(core_genes),
        " | heterogeneity_flagged_core=",
        sum(out$heterogeneity_flag, na.rm = TRUE))

message("\nDONE: Step 06B (anchor heterogeneity) completed.")