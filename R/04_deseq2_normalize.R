suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 04_deseq2_normalize.R <GSE_ID> [project_root]\n",
       "Example: Rscript 04_deseq2_normalize.R GSE264344 /orange/qsong1/Yansheng")
}

dataset_id   <- args[1]
project_root <- ifelse(length(args) >= 2, args[2], getwd())

counts_path <- file.path(project_root, "03_counts", dataset_id,
                         "featurecounts", "validation", "gene_counts_clean.tsv")
design_path <- file.path(project_root, "00_metadata", dataset_id,
                         paste0(dataset_id, "_design.tsv"))

out_dir <- file.path(project_root, "04_de", dataset_id, "normalized")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)
if (!file.exists(design_path)) stop("Design file not found: ", design_path)

message("Dataset: ", dataset_id)
message("Project root: ", project_root)
message("Counts: ", counts_path)
message("Design: ", design_path)

# -------------------------
# Load counts (genes x samples)
# -------------------------
counts <- read_tsv(counts_path, show_col_types = FALSE)

message("Counts columns: ", paste(colnames(counts), collapse = ", "))

# robust gene id column detection
gene_col <- dplyr::case_when(
  "gene_id" %in% names(counts) ~ "gene_id",
  "Geneid"  %in% names(counts) ~ "Geneid",
  "GeneID"  %in% names(counts) ~ "GeneID",
  "gene"    %in% names(counts) ~ "gene",
  TRUE                         ~ names(counts)[1]
)

message("Using gene id column: ", gene_col)

counts <- counts %>% dplyr::rename(gene_id = !!gene_col)

count_mat <- counts %>%
  column_to_rownames("gene_id") %>%
  as.data.frame()

# force numeric matrix (featureCounts outputs are typically integer-like)
count_mat[] <- lapply(count_mat, function(x) as.numeric(x))
count_mat <- as.matrix(count_mat)

if (anyNA(count_mat)) {
  stop("NA detected in count matrix after numeric coercion. Check gene_counts_clean.tsv for non-numeric entries.")
}

# DESeq2 expects integer counts
storage.mode(count_mat) <- "integer"

# -------------------------
# Load design (must include sample_id)
# -------------------------
design <- read_tsv(design_path, show_col_types = FALSE)

if (!("sample_id" %in% colnames(design))) {
  stop("design file must contain a 'sample_id' column.")
}

design <- design %>%
  filter(sample_id %in% colnames(count_mat)) %>%
  arrange(match(sample_id, colnames(count_mat)))

if (nrow(design) == 0) stop("No sample_id in design matched count matrix columns.")

if (!all(design$sample_id == colnames(count_mat))) {
  stop("Sample order mismatch: design$sample_id must match count matrix columns.")
}

design <- design %>% column_to_rownames("sample_id")

# -------------------------
# DESeq2: size-factor normalization only (no inference)
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = design,
  design    = ~ 1
)

# optional: drop very low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- estimateSizeFactors(dds)

# -------------------------
# Export normalized counts + size factors
# -------------------------
norm_counts <- counts(dds, normalized = TRUE)

norm_df <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene_id")

write_tsv(norm_df, file.path(out_dir, "gene_counts_normalized.tsv"))

sf_df <- tibble(
  sample_id   = colnames(dds),
  size_factor = sizeFactors(dds)
)

write_tsv(sf_df, file.path(out_dir, "size_factors.tsv"))

message("Done.")
message("Wrote: ", file.path(out_dir, "gene_counts_normalized.tsv"))
message("Wrote: ", file.path(out_dir, "size_factors.tsv"))
