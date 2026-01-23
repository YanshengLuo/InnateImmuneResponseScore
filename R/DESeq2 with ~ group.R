library(DESeq2)
library(readr)
library(dplyr)

#read_inputs
dataset_id <- "GSE264344"
tissue <- "Blood"
project_root <- "C:/Users/john/Desktop/IMRS_Project"

# counts: genes x samples
counts <- read_tsv(
  file.path(project_root, "03_counts",
            paste0(dataset_id, "_", tissue, "_counts.tsv")),
  show_col_types = FALSE
)

# design matrix
design <- read_tsv(
  file.path(project_root, "00_metadata",
            paste0(dataset_id, "_", tissue, "_design.tsv")),
  show_col_types = FALSE
)

#prepare DESeq2 objects
# convert counts to matrix
count_mat <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# ensure order matches
design <- design %>%
  filter(sample_id %in% colnames(count_mat)) %>%
  arrange(match(sample_id, colnames(count_mat)))

stopifnot(all(design$sample_id == colnames(count_mat)))

#create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = design,
  design = ~ group
)

# baseline reference
dds$group <- relevel(dds$group, ref = "baseline_0")

dds <- DESeq(dds)


##Extract contrasts **
groups <- levels(dds$group)

# delivery groups only
contrast_groups <- groups[groups != "baseline_0"]

out_dir <- file.path(project_root, "04_de", dataset_id, tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (g in contrast_groups) {
  
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
