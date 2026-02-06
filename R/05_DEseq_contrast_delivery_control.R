#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 5: Within-dataset DE contrast extraction (delivery vs control)
#
# KEEP THIS IDEA (for any AI / future edits):
# - IMRS is CONTRAST-driven, not pooled-expression-driven.
# - We NEVER compare tissue vs tissue as a biological contrast.
# - For each dataset, estimate DELIVERY vs CONTROL within the SAME tissue and timepoint.
# - Tissue/time are metadata labels used for subsetting, not the main effect.
# - Time bins:
#     <24h   = "anchor" contrasts
#     24-72h = "calibration" contrasts
# - Output unit: ONE DE table per (dataset, tissue, time_bin+time, delivery_vs_control)
#
# Inputs (your layout):
# - Normalized counts:
#   /orange/qsong1/Yansheng/04_de/GSE*/normalized/gene_counts_normalized.tsv
#   (also checks typo: gene_counts_normailzed.tsv)
# - Size factors:
#   /orange/qsong1/Yansheng/04_de/GSE*/normalized/size_factors.tsv
# - Design:
#   /orange/qsong1/Yansheng/00_metadata/GSE*/GSE*_design.tsv
#
# IMPORTANT:
# - DESeq2 expects RAW integer counts.
# - You currently have normalized counts + size factors.
# - Since normalized = raw / sizeFactor, we reconstruct approx raw:
#       raw ≈ normalized * sizeFactor
#   then round to integers.
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# -------------------------
# Args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript 05_deseq2_contrasts.R <GSE_ID> [project_root]\n",
    "Example: Rscript 05_deseq2_contrasts.R GSE264344 /orange/qsong1/Yansheng\n"
  )
}
dataset_id   <- args[1]
project_root <- ifelse(length(args) >= 2, args[2], "/orange/qsong1/Yansheng")

# -------------------------
# Paths
# -------------------------
norm_dir <- file.path(project_root, "04_de", dataset_id, "normalized")

counts_candidates <- c(
  file.path(norm_dir, "gene_counts_normalized.tsv"),
  file.path(norm_dir, "gene_counts_normailzed.tsv") # handle typo
)
counts_path <- counts_candidates[file.exists(counts_candidates)][1]
if (is.na(counts_path)) stop("Normalized counts file not found in: ", norm_dir)

sf_path <- file.path(norm_dir, "size_factors.tsv")
if (!file.exists(sf_path)) stop("Size factors file not found: ", sf_path)

design_candidates <- c(
  file.path(project_root, "00_metadata", dataset_id, paste0(dataset_id, "_design.tsv")),
  file.path(project_root, "00_metadata", paste0(dataset_id, "_design.tsv"))
)
design_path <- design_candidates[file.exists(design_candidates)][1]
if (is.na(design_path)) stop("Design file not found for dataset: ", dataset_id)

out_root <- file.path(project_root, "04_de", dataset_id, "deseq2_contrasts")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

message("Dataset: ", dataset_id)
message("Project root: ", project_root)
message("Normalized counts: ", counts_path)
message("Size factors: ", sf_path)
message("Design: ", design_path)
message("Output root: ", out_root)

# ----------------------------
# Column synonyms (IMRS standardized dictionary)
# ----------------------------
syn <- list(
  srr_id        = c("Run", "run", "SRR", "srr", "srr_id", "srr_accession"),
  biosample_id  = c("BioSample", "biosample", "BioSample_ID", "biosample_id"),
  experiment_id = c("Experiment", "experiment", "SRX", "srx", "experiment_id"),
  tissue        = c("tissue", "organ", "Organ", "tissue_type"),
  source_name   = c("source_name", "source_name_ch1", "SourceName", "source"),
  timepoint     = c("timepoint", "time_point", "time", "Time",
                    "collection_time", "hours_post", "hrs_post"),
  delivery_type = c("vaccine_type", "vaccine", "treatment", "treat",
                    "platform", "vector", "delivery", "group", "arm")
)

# Robust column picker
pick_col <- function(df, choices) {
  nm <- names(df)
  hit <- choices[choices %in% nm]
  if (length(hit) >= 1) return(hit[1])

  # case-insensitive fallback
  nm_low <- tolower(nm)
  choices_low <- tolower(choices)
  idx <- match(choices_low, nm_low)
  idx <- idx[!is.na(idx)][1]
  if (!is.na(idx)) return(nm[idx])

  return(NA_character_)
}

# -------------------------
# Load design
# -------------------------
design <- read_tsv(design_path, show_col_types = FALSE)

if (!("sample_id" %in% names(design))) {
  stop("design file must contain a 'sample_id' column.")
}

col_map <- list(
  srr_id        = pick_col(design, syn$srr_id),
  biosample_id  = pick_col(design, syn$biosample_id),
  experiment_id = pick_col(design, syn$experiment_id),
  tissue        = pick_col(design, syn$tissue),
  source_name   = pick_col(design, syn$source_name),
  timepoint     = pick_col(design, syn$timepoint),
  delivery_type = pick_col(design, syn$delivery_type)
)

message("Detected columns (col_map):")
print(col_map)

tissue_col <- col_map$tissue
time_col   <- col_map$timepoint
cond_col   <- col_map$delivery_type
batch_col  <- pick_col(design, c("batch", "Batch", "run", "lane")) # optional (not always used)

if (is.na(cond_col)) {
  stop("Could not detect delivery/condition column using dictionary. ",
       "Add a column like 'group'/'treatment'/'delivery' or extend syn$delivery_type.")
}

# -------------------------
# Load normalized counts
# -------------------------
norm_df <- read_tsv(counts_path, show_col_types = FALSE)

# Robust gene id detection
gene_col <- dplyr::case_when(
  "gene_id" %in% names(norm_df) ~ "gene_id",
  "Geneid"  %in% names(norm_df) ~ "Geneid",
  "GeneID"  %in% names(norm_df) ~ "GeneID",
  "gene"    %in% names(norm_df) ~ "gene",
  TRUE                          ~ names(norm_df)[1]
)
norm_df <- norm_df %>% rename(gene_id = !!gene_col)

norm_mat <- norm_df %>%
  column_to_rownames("gene_id") %>%
  as.data.frame()

message("Normalized matrix: ", nrow(norm_mat), " genes x ", ncol(norm_mat), " samples")

# -------------------------
# Load size factors
# -------------------------
sf_df <- read_tsv(sf_path, show_col_types = FALSE)

# Accept either (sample_id,size_factor) or (sample_id,sizeFactors) etc.
sf_col <- pick_col(sf_df, c("size_factor", "sizeFactor", "size_factors", "sizeFactors"))
if (!("sample_id" %in% names(sf_df)) || is.na(sf_col)) {
  stop("size_factors.tsv must contain 'sample_id' and a size factor column (e.g., size_factor).")
}
sf_df <- sf_df %>% rename(size_factor = !!sf_col)

# -------------------------
# Restrict to shared samples and align order
# -------------------------
shared_samples <- Reduce(intersect, list(colnames(norm_mat), design$sample_id, sf_df$sample_id))
if (length(shared_samples) < 4) {
  stop("Not enough shared samples among normalized counts, design, and size_factors (need >=4).")
}

norm_mat <- norm_mat[, shared_samples, drop = FALSE]
design2  <- design %>%
  filter(sample_id %in% shared_samples) %>%
  arrange(match(sample_id, shared_samples))
sf_df2   <- sf_df %>%
  filter(sample_id %in% shared_samples) %>%
  arrange(match(sample_id, shared_samples))

stopifnot(all(design2$sample_id == shared_samples))
stopifnot(all(sf_df2$sample_id == shared_samples))

# -------------------------
# Reconstruct approx raw counts: raw ≈ normalized * size_factor
# -------------------------
norm_mat[] <- lapply(norm_mat, as.numeric)
norm_mat <- as.matrix(norm_mat)
if (anyNA(norm_mat)) stop("NA in normalized matrix after coercion. Check normalized TSV.")

sf <- sf_df2$size_factor
names(sf) <- sf_df2$sample_id

raw_mat <- sweep(norm_mat, 2, sf[colnames(norm_mat)], `*`)
raw_mat <- round(raw_mat)
raw_mat[raw_mat < 0] <- 0
storage.mode(raw_mat) <- "integer"

# -------------------------
# Build standard metadata columns
# -------------------------
design2 <- design2 %>%
  mutate(
    .delivery_raw = as.character(.data[[cond_col]])
  )

if (!is.na(tissue_col)) {
  design2 <- design2 %>% mutate(.tissue = as.character(.data[[tissue_col]]))
} else {
  design2 <- design2 %>% mutate(.tissue = "NA_TISSUE")
}

# Parse timepoint hours (best-effort)
if (!is.na(time_col)) {
  time_raw <- as.character(design2[[time_col]])
  time_num <- suppressWarnings(as.numeric(time_raw))
  if (all(is.na(time_num))) {
    # extract numeric from strings like "6h", "24 hr", "48hours"
    time_num <- suppressWarnings(as.numeric(str_extract(time_raw, "\\d+\\.?\\d*")))
  }
  design2 <- design2 %>% mutate(.time_h = time_num)
} else {
  design2 <- design2 %>% mutate(.time_h = NA_real_)
}

design2 <- design2 %>%
  mutate(
    .time_bin = case_when(
      !is.na(.time_h) & .time_h < 24                    ~ "anchor_lt24h",
      !is.na(.time_h) & .time_h >= 24 & .time_h <= 72   ~ "calib_24to72h",
      TRUE                                              ~ "time_unknown_or_outside"
    ),
    .time_label = ifelse(!is.na(.time_h), paste0(.time_h, "h"), "NAh")
  )

# -------------------------
# Identify controls by keyword (EDIT HERE if your labels differ)
# -------------------------
# This checks the delivery/group string for control-like labels.
control_regex <- "(^|[_ -])(control|ctrl|pbs|vehicle|mock|naive|saline|untreated|sham|empty)([_ -]|$)"

design2 <- design2 %>%
  mutate(
    .is_control = str_detect(str_to_lower(.delivery_raw), control_regex),
    .condition_clean = case_when(
      .is_control ~ "CONTROL",
      TRUE        ~ paste0("DELIVERY_", make.names(.delivery_raw))
    )
  )

message("Control detection summary:")
print(table(design2$.condition_clean, useNA = "ifany"))

# -------------------------
# Contrast extraction settings
# -------------------------
min_n_per_group <- 2  # require >=2 control and >=2 delivery per contrast
gene_min_total  <- 10 # filter genes with total counts < 10 in subset

# -------------------------
# Function: run a single delivery vs control contrast for a subset
# -------------------------
run_contrast <- function(sub_samples, tissue, time_bin, time_label, delivery_level) {
  sub_design <- design2 %>%
    filter(sample_id %in% sub_samples) %>%
    mutate(
      condition_simple = factor(.condition_clean, levels = c("CONTROL", delivery_level))
    )

  tab <- table(sub_design$condition_simple)
  if (!all(c("CONTROL", delivery_level) %in% names(tab))) return(NULL)
  if (any(tab[c("CONTROL", delivery_level)] < min_n_per_group)) return(NULL)

  sub_counts <- raw_mat[, sub_design$sample_id, drop = FALSE]

  coldata <- sub_design %>%
    select(sample_id, condition_simple, .tissue, .time_h, .time_bin, .time_label, .delivery_raw) %>%
    column_to_rownames("sample_id")

  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData   = coldata,
    design    = ~ condition_simple
  )

  dds <- dds[rowSums(counts(dds)) >= gene_min_total, ]
  if (nrow(dds) < 50) {
    # too few genes left is usually a sign something is wrong in subset
    return(NULL)
  }

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition_simple", delivery_level, "CONTROL"))

  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      dataset_id = dataset_id,
      tissue     = tissue,
      time_bin   = time_bin,
      time_label = time_label,
      contrast   = paste0(delivery_level, "_vs_CONTROL"),
      n_control  = as.integer(tab["CONTROL"]),
      n_delivery = as.integer(tab[delivery_level])
    )

  safe_tissue <- make.names(tissue)
  safe_time   <- make.names(paste(time_bin, time_label, sep="__"))
  safe_deliv  <- make.names(delivery_level)

  out_dir <- file.path(out_root, safe_tissue, safe_time)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(
    out_dir,
    paste0(dataset_id, "__", safe_tissue, "__", safe_time, "__", safe_deliv, "_vs_CONTROL.tsv")
  )

  write_tsv(res_df, out_file)
  message("Wrote: ", out_file)
  out_file
}

# -------------------------
# Main loop: tissue x time -> delivery vs control
# -------------------------
all_tissues <- unique(design2$.tissue)
all_times   <- unique(paste(design2$.time_bin, design2$.time_label, sep="__"))

written <- character(0)

for (t in all_tissues) {
  for (tt in all_times) {
    parts <- str_split(tt, "__", simplify = TRUE)
    time_bin   <- parts[1]
    time_label <- parts[2]

    sub <- design2 %>%
      filter(.tissue == t) %>%
      filter(paste(.time_bin, .time_label, sep="__") == tt)

    if (nrow(sub) == 0) next
    if (!any(sub$.condition_clean == "CONTROL")) next

    delivery_levels <- setdiff(unique(sub$.condition_clean), "CONTROL")
    if (length(delivery_levels) == 0) next

    for (dl in delivery_levels) {
      sub_samples <- sub %>%
        filter(.condition_clean %in% c("CONTROL", dl)) %>%
        pull(sample_id)

      out <- run_contrast(sub_samples, t, time_bin, time_label, dl)
      if (!is.null(out)) written <- c(written, out)
    }
  }
}

# -------------------------
# Write an index of outputs
# -------------------------
index_file <- file.path(out_root, paste0(dataset_id, "__contrast_index.tsv"))
if (length(written) > 0) {
  tibble(result_file = written) %>% write_tsv(index_file)
  message("Done. Index: ", index_file)
} else {
  message("Done. No contrasts were written.")
  message("Most common reasons:")
  message("  - control labels not matched by control_regex (edit control_regex)")
  message("  - timepoint column missing/unparseable (time becomes NAh; still runs but may collapse)")
  message("  - not enough replicates per group (min_n_per_group = ", min_n_per_group, ")")
}
