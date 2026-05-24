#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(biomaRt)
})

# =========================
# USER SETTINGS
# =========================
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

weights_file <- file.path(project_root, "05_score", "anchors", "gene_weights.tsv")
out_dir <- file.path(project_root, "05_score", "human_transfer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ortholog_file <- file.path(out_dir, "mouse_human_ortholog_table.tsv")
human_weights_file <- file.path(out_dir, "human_gene_weights.tsv")
coverage_file <- file.path(out_dir, "ortholog_coverage.tsv")

# =========================
# HELPERS
# =========================
strip_version <- function(x) {
  sub("\\..*$", "", x)
}

# =========================
# READ WEIGHTS
# =========================
w <- read_tsv(weights_file, show_col_types = FALSE)

if (!all(c("gene", "weight") %in% colnames(w))) {
  stop("gene_weights.tsv must contain at least: gene, weight")
}

w <- w %>%
  mutate(
    mouse_gene_id_version = gene,
    mouse_gene_id = strip_version(gene)
  ) %>%
  distinct(mouse_gene_id, .keep_all = TRUE)

n_mouse_total <- nrow(w)

# =========================
# CONNECT TO ENSEMBL
# =========================
mouse_mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

# =========================
# QUERY MOUSE -> HUMAN ORTHOLOGS
# =========================
# We pull:
# - mouse Ensembl gene ID
# - human ortholog Ensembl gene ID
# - orthology type
# - optional symbols for readability
ortho_raw <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type"
  ),
  filters = "ensembl_gene_id",
  values = unique(w$mouse_gene_id),
  mart = mouse_mart
)

# =========================
# STRICT 1:1 FILTER
# =========================
ortho_1to1 <- ortho_raw %>%
  rename(
    mouse_gene_id = ensembl_gene_id,
    mouse_symbol = external_gene_name,
    human_gene_id = hsapiens_homolog_ensembl_gene,
    human_symbol = hsapiens_homolog_associated_gene_name,
    orthology_type = hsapiens_homolog_orthology_type
  ) %>%
  filter(
    !is.na(human_gene_id),
    human_gene_id != "",
    orthology_type == "ortholog_one2one"
  ) %>%
  distinct(mouse_gene_id, .keep_all = TRUE) %>%
  distinct(human_gene_id, .keep_all = TRUE)

# =========================
# JOIN WEIGHTS
# =========================
human_weights <- w %>%
  inner_join(ortho_1to1, by = "mouse_gene_id") %>%
  transmute(
    mouse_gene_id_version,
    mouse_gene_id,
    mouse_symbol,
    human_gene_id,
    human_symbol,
    orthology_type,
    weight
  ) %>%
  arrange(desc(abs(weight)))

n_mapped <- nrow(human_weights)
coverage <- n_mapped / n_mouse_total

coverage_tbl <- tibble(
  mouse_weight_genes_total = n_mouse_total,
  one_to_one_mapped_genes = n_mapped,
  ortholog_coverage = coverage
)

# =========================
# WRITE OUTPUTS
# =========================
write_tsv(ortho_1to1, ortholog_file)
write_tsv(human_weights, human_weights_file)
write_tsv(coverage_tbl, coverage_file)

message("Wrote: ", ortholog_file)
message("Wrote: ", human_weights_file)
message("Wrote: ", coverage_file)
message("Ortholog coverage = ", round(coverage, 4))