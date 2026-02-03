library(DESeq2)
library(readr)
library(dplyr)
library(tibble)  # needed for column_to_rownames / rownames_to_column

# ============================================================
# 05_DESeq2_DE.R
#
# PURPOSE (IMRS workflow mapping)
#   INPUT:
#     - gene count matrix (genes x samples), raw integer counts
#     - design matrix defining groups (baseline vs delivery conditions)
#
#   OUTPUT:
#     - One DE table per delivery group vs baseline_0:
#         04_de/<dataset_id>/<tissue>/DE_<group>_vs_baseline_0.tsv
#
# IMRS WORKFLOW STEP MAPPING
#   Step 4 (Per-dataset normalization):
#     - DESeq2 size-factor normalization occurs automatically inside DESeq(dds)
#     - This makes samples comparable within THIS dataset only
#
#   Step 5 (Differential expression):
#     - Fit negative binomial GLM and test: delivery - baseline_0
#
# IMPORTANT RULES (IMRS)
#   - Do NOT merge counts across datasets here
#   - Do NOT batch-correct across studies here
#   - Each dataset+tissue analyzed independently; cross-dataset logic happens later
# ============================================================

# -------------------------
# read_inputs
# -------------------------
dataset_id <- "GSE264344"
tissue <- "Blood"
project_root <- "C:/Users/john/Desktop/IMRS_Project"

# counts: genes x samples
# expectation: raw integer counts (NOT TPM/CPM/normalized)
counts <- read_tsv(
  file.path(project_root, "03_counts",
            paste0(dataset_id, "_", tissue, "_counts.tsv")),
  show_col_types = FALSE
)

# design matrix: sample metadata defining the model groups
# expectation: contains at least:
#   - sample_id (matching count matrix column names)
#   - group (baseline_0 + one or more delivery groups)
design <- read_tsv(
  file.path(project_root, "00_metadata",
            paste0(dataset_id, "_", tissue, "_design.tsv")),
  show_col_types = FALSE
)

# -------------------------
# prepare DESeq2 objects
# -------------------------

# convert counts table -> matrix with gene_id as rownames
# result: count_mat[gene_id, sample_id]
count_mat <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# -------------------------
# align design to count matrix columns
# -------------------------
# This ensures design rows are in the same order as count_mat columns.
# If misaligned, DESeq2 would associate the wrong metadata with samples.
design <- design %>%
  filter(sample_id %in% colnames(count_mat)) %>%
  arrange(match(sample_id, colnames(count_mat)))

stopifnot(all(design$sample_id == colnames(count_mat)))

# Make the DESeq2 mapping explicit:
# DESeq2 uses rownames(colData) to match samples to countData columns.
design_df <- as.data.frame(design)
rownames(design_df) <- design_df$sample_id

# -------------------------
# create DESeq2 dataset
# -------------------------
# Model: expression ~ group
# Interpretation: we estimate group-specific shifts relative to the reference.
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = design_df,
  design = ~ group
)

# -------------------------
# baseline reference level
# -------------------------
# Setting baseline_0 as reference makes log2FoldChange interpretable as:
#   log2FC > 0  => higher expression in delivery group than baseline_0
#   log2FC < 0  => lower expression in delivery group than baseline_0
dds$group <- factor(dds$group)
dds$group <- relevel(dds$group, ref = "baseline_0")

# ============================================================
# STEP 4 (IMRS): Per-dataset normalization (size factors)
# ============================================================
# When you run DESeq(dds), DESeq2 performs (in order):
#   1) estimateSizeFactors()  <-- IMRS Step 4 normalization
#   2) estimateDispersions()
#   3) fit GLM + Wald tests   <-- IMRS Step 5 differential expression
#
# Optional (if you want an explicit Step 4 output file):
#   dds_sf <- estimateSizeFactors(dds)
#   norm_counts <- counts(dds_sf, normalized = TRUE)
#   write_tsv(as.data.frame(norm_counts) %>% rownames_to_column("gene_id"),
#             file.path(project_root, "03_counts",
#                       paste0(dataset_id, "_", tissue, "_normcounts.tsv")))
# ============================================================

dds <- DESeq(dds)

# -------------------------
# extract contrasts
# -------------------------
groups <- levels(dds$group)

# delivery groups only: compare each non-baseline group against baseline_0
contrast_groups <- groups[groups != "baseline_0"]

out_dir <- file.path(project_root, "04_de", dataset_id, tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (g in contrast_groups) {

  # DESeq2 contrast definition:
  #   contrast = c("group", g, "baseline_0")
  #
  # Meaning:
  #   log2FoldChange = log2( mean expression in group g / mean expression in baseline_0 )
  #
  # Therefore:
  #   log2FoldChange > 0  => gene is UP in delivery group g vs baseline_0
  #   log2FoldChange < 0  => gene is DOWN in delivery group g vs baseline_0
  #
  # pvalue / padj test whether the group difference is non-zero (Wald test).
  res <- results(
    dds,
    contrast = c("group", g, "baseline_0"),
    alpha = 0.05
  )

  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      dataset_id = dataset_id,
      tissue = tissue,
      contrast = paste0(g, "_vs_baseline_0")
    ) %>%
    select(
      gene_id,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      dataset_id,
      tissue,
      contrast
    )

  write_tsv(
    res_tbl,
    file.path(out_dir, paste0("DE_", g, "_vs_baseline_0.tsv"))
  )
}

# ============================================================
# WHAT TO DO NEXT (IMRS workflow continuation)
# ============================================================
# After this script finishes for ALL datasets/tissues:
#
# Next Step (IMRS Step 6): Anchor-based reproducibility filtering
#   - Collect DE tables from anchor datasets only (mouse LNP/LV/AdV anchors)
#   - Build the Core Gene Set using majority-consensus rules:
#       * gene upregulated in >= 2 anchors
#       * consistent direction among those anchors
#       * effect-size threshold (e.g., |log2FC| >= 0.5)
#       * statistical support in >= 1 anchor (padj <= threshold)
#
# Output of Step 6:
#   model/core_gene_set.tsv
#
# Then (IMRS Step 7): Gene weighting across anchors
#   - Use meta-analytic aggregation of log2FC (optionally inverse-variance)
#   - Cap extreme weights and down-weight high-variance genes
#
# Then (IMRS Steps 8–10): scoring per dataset (z-scoring + weighted sum)
# ============================================================
