# ============================================================
# IMRS Step 5 (Windows, interactive):
# DESeq2 contrasts from VERIFIED split design folders
#
# Split designs:
#   C:/Users/john/Desktop/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/
#     verified_metadata/splited/<DATASET>_design/*.tsv
#
# Normalized inputs:
#   C:/Users/john/Desktop/IMRS_Project/04_de/<DATASET>/normalized/gene_counts_normalized.tsv
#   C:/Users/john/Desktop/IMRS_Project/04_de/<DATASET>/normalized/size_factors.tsv
#
# Outputs:
#   C:/Users/john/Desktop/IMRS_Project/04_de/comparison/<DATASET>/deseq2_contrasts/*.tsv
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

project_root <- "C:/Users/john/Desktop/IMRS_Project"
de_root <- file.path(project_root, "04_de")

split_design_root <- file.path(
  project_root,
  "Hypergator_scripts", "InnateImmuneResponseScore",
  "verified_metadata", "splited"
)

out_base <- file.path(project_root, "04_de", "comparison")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

min_n_per_group <- 2
gene_min_total  <- 10

# -------- helpers --------
is_real_split_design <- function(path) {
  bn <- tolower(basename(path))
  if (!str_detect(bn, "\\.tsv$")) return(FALSE)
  if (str_detect(bn, "contrast_index|group_map")) return(FALSE)
  str_detect(bn, "design__t=")
}

pick_first <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit) > 0) return(hit[1])
  nm_low <- tolower(nm); cand_low <- tolower(candidates)
  idx <- match(cand_low, nm_low)
  idx <- idx[!is.na(idx)][1]
  if (!is.na(idx)) return(nm[idx])
  NA_character_
}

read_table_robust <- function(path) {
  # tab
  df <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  # comma
  df <- tryCatch(read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  # whitespace
  df <- tryCatch(read_table(path, show_col_types = FALSE, progress = FALSE),
                 error = function(e) NULL)
  if (!is.null(df) && ncol(df) > 1) return(df)

  stop("Failed to read (delimiter/header issue): ", path)
}

# This is the key fix: build a numeric matrix with explicit colnames.
make_counts_matrix <- function(df) {
  gene_col <- dplyr::case_when(
    "gene_id" %in% names(df) ~ "gene_id",
    "Geneid"  %in% names(df) ~ "Geneid",
    "GeneID"  %in% names(df) ~ "GeneID",
    "gene"    %in% names(df) ~ "gene",
    TRUE                     ~ names(df)[1]
  )
  df <- df %>% rename(gene_id = !!gene_col)

  if (ncol(df) < 3) {
    stop("Counts table has <3 columns (need gene_id + >=2 samples). Columns: ",
         paste(names(df), collapse = ", "))
  }

  sample_cols <- setdiff(names(df), "gene_id")
  gene_ids <- as.character(df$gene_id)

  mat <- as.matrix(df[, sample_cols])
  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- gene_ids
  colnames(mat) <- trimws(sample_cols)

  mat
}

# -------- per dataset runner --------
run_one_dataset <- function(dataset_id) {

  message("\n==============================")
  message("DATASET: ", dataset_id)
  message("==============================")

  # normalized
  norm_dir <- file.path(de_root, dataset_id, "normalized")
  counts_candidates <- c(
    file.path(norm_dir, "gene_counts_normalized.tsv"),
    file.path(norm_dir, "gene_counts_normailzed.tsv")
  )
  counts_path <- counts_candidates[file.exists(counts_candidates)][1]
  sf_path <- file.path(norm_dir, "size_factors.tsv")

  if (is.na(counts_path) || !file.exists(sf_path)) {
    message("Skip: missing counts or size_factors in ", norm_dir)
    return(invisible(NULL))
  }

  # designs
  design_dir <- file.path(split_design_root, paste0(dataset_id, "_design"))
  if (!dir.exists(design_dir)) {
    message("Skip: missing split design folder: ", design_dir)
    return(invisible(NULL))
  }

  design_files <- list.files(design_dir, pattern = "\\.tsv$", full.names = TRUE)
  design_files <- design_files[sapply(design_files, is_real_split_design)]
  if (length(design_files) == 0) {
    message("Skip: no per-contrast design TSVs found in: ", design_dir)
    return(invisible(NULL))
  }

  message("Found split design files: ", length(design_files))
  message("DEBUG counts_path = ", counts_path)
  message("DEBUG sf_path     = ", sf_path)

  # ---- load counts (robust) and build matrix (FIXED) ----
  norm_df <- read_table_robust(counts_path)
  norm_mat <- make_counts_matrix(norm_df)

  # ---- load size factors ----
  sf_raw <- read_table_robust(sf_path)
  id_col <- pick_first(names(sf_raw), c("sample_id","sample","Run","run","SRR","srr","id"))
  sf_col <- pick_first(names(sf_raw), c("size_factor","sizeFactor","sizeFactors","size_factors"))

  if (is.na(id_col) || is.na(sf_col)) {
    stop("size_factors.tsv missing id/size_factor columns: ", sf_path,
         "\nColumns: ", paste(names(sf_raw), collapse = ", "))
  }

  sf_df <- sf_raw %>%
    transmute(
      sample_id = trimws(as.character(.data[[id_col]])),
      size_factor = suppressWarnings(as.numeric(.data[[sf_col]]))
    )

  # ---- shared samples ----
  message("DEBUG counts head = ", paste(head(colnames(norm_mat), 8), collapse = ", "))
  message("DEBUG sf head     = ", paste(head(sf_df$sample_id, 8), collapse = ", "))

  shared <- intersect(colnames(norm_mat), sf_df$sample_id)
  message("DEBUG shared n    = ", length(shared))

  if (length(shared) < 4) {
    message("Skip: too few shared samples between counts and size_factors for ", dataset_id)
    return(invisible(NULL))
  }

  norm_mat <- norm_mat[, shared, drop = FALSE]
  sf_df <- sf_df %>% filter(sample_id %in% shared) %>% arrange(match(sample_id, shared))

  # reconstruct raw
  sf_vec <- sf_df$size_factor
  names(sf_vec) <- sf_df$sample_id

  raw_mat <- sweep(norm_mat, 2, sf_vec[colnames(norm_mat)], `*`)
  raw_mat <- round(raw_mat)
  raw_mat[raw_mat < 0] <- 0
  storage.mode(raw_mat) <- "integer"

  # outputs
  out_root <- file.path(out_base, dataset_id, "deseq2_contrasts")
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  written <- character(0)

  for (df_path in design_files) {

    d <- read_table_robust(df_path)

    if (!all(c("sample_id","condition_simple") %in% names(d))) {
      message("Skip design (missing sample_id/condition_simple): ", basename(df_path))
      next
    }

    d <- d %>%
      mutate(
        sample_id = trimws(as.character(sample_id)),
        condition_simple = as.character(condition_simple)
      )

    if (!all(c("CONTROL","DELIVERY") %in% unique(d$condition_simple))) {
      message("Skip design (missing CONTROL/DELIVERY): ", basename(df_path))
      next
    }

    tab <- table(d$condition_simple)
    if (any(tab[c("CONTROL","DELIVERY")] < min_n_per_group)) {
      message("Skip design (<", min_n_per_group, " per group): ", basename(df_path))
      next
    }

    samples <- intersect(d$sample_id, colnames(raw_mat))
    if (length(samples) < 4) {
      message("Skip design (too few samples found in counts): ", basename(df_path))
      next
    }

    d2 <- d %>%
      filter(sample_id %in% samples) %>%
      mutate(condition_simple = factor(condition_simple, levels = c("CONTROL","DELIVERY"))) %>%
      arrange(match(sample_id, samples))

    counts_sub <- raw_mat[, d2$sample_id, drop = FALSE]

    coldata <- d2 %>%
      select(sample_id, condition_simple) %>%
      column_to_rownames("sample_id")

    dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ condition_simple)
    dds <- dds[rowSums(counts(dds)) >= gene_min_total, ]
    if (nrow(dds) < 50) {
      message("Skip design (too few genes after filtering): ", basename(df_path))
      next
    }

    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition_simple","DELIVERY","CONTROL"))

    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      mutate(
        dataset_id  = dataset_id,
        design_file = basename(df_path),
        n_control   = as.integer(tab["CONTROL"]),
        n_delivery  = as.integer(tab["DELIVERY"])
      )

    out_file <- file.path(
      out_root,
      paste0(dataset_id, "__", tools::file_path_sans_ext(basename(df_path)), "__DE.tsv")
    )

    write_tsv(res_df, out_file)
    written <- c(written, out_file)
    message("Wrote: ", out_file)
  }

  index_file <- file.path(out_root, paste0(dataset_id, "__contrast_index.tsv"))
  if (length(written) > 0) {
    tibble(result_file = written) %>% write_tsv(index_file)
    message("Index: ", index_file)
  } else {
    message("No contrasts written for ", dataset_id)
  }

  invisible(NULL)
}

# -------- MAIN LOOP --------
dataset_ids <- list.dirs(de_root, recursive = FALSE, full.names = FALSE)
dataset_ids <- dataset_ids[grepl("^GSE", dataset_ids)]
print(dataset_ids)

for (ds in dataset_ids) run_one_dataset(ds)

message("\nAll DE contrasts finished. Outputs under: ", out_base)
