# ============================================================
# IMRS Step 06A — Core gene sets (LOCKED mouse datasets)
# - Builds core gene set separately for:
#     (1) ANCHOR phase      (H <= 24)   -> Step 5 output folder: .../anchor/
#     (2) CALIBRATION phase (24 < H <72)-> Step 5 output folder: .../calibration/
#
# Reads Step 5 workflow DE tables:
#   gene, log2FC, SE, stat, pval, FDR
#
# Outputs:
#   <project_root>/05_score/anchors/core_gene_set.tsv
#   <project_root>/05_score/calibration/core_gene_set.tsv
#   plus small contrast-count audits per phase
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

# Locked mouse datasets (pre-registered)
LOCKED_DATASETS_MOUSE <- c("GSE39129", "GSE167521", "GSE264344")

# Thresholds
fdr_cutoff <- 0.05
lfc_cutoff <- 1.0
within_dataset_support_cutoff <- 0.50   # gene passes in >=50% contrasts within dataset
consensus_k_fraction          <- 2/3    # gene supported by >= ceil(K_present*2/3) datasets

# Which phases to build
PHASES <- c("anchor", "calibration")

# -------------------------
# Helpers
# -------------------------
read_de_workflow <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("gene", "log2FC", "SE", "stat", "pval", "FDR")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in: ", path, "\nMissing: ", paste(missing, collapse = ", "))
  }
  df %>%
    transmute(
      gene_id = as.character(gene),
      log2FC  = suppressWarnings(as.numeric(log2FC)),
      FDR     = suppressWarnings(as.numeric(FDR))
    )
}

get_dataset_from_path <- function(path) {
  # .../comparison/<DATASET>/deseq2_contrasts/<phase>/<file>
  parts <- strsplit(normalizePath(path, winslash = "/"), "/")[[1]]
  i <- which(parts == "comparison")
  if (length(i) == 1 && length(parts) >= i + 1) return(parts[i + 1])
  NA_character_
}

# Build core gene set for one phase folder (anchor or calibration)
build_core_for_phase <- function(phase) {
  stopifnot(phase %in% c("anchor", "calibration"))

  out_root <- file.path(project_root, "05_score", phase)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  # Collect workflow files from this phase
  all_files <- list.files(
    de_comp_root,
    pattern = "__DE_workflow\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  phase_files <- all_files[str_detect(all_files, paste0("deseq2_contrasts[\\\\/]+", phase, "[\\\\/]"))]

  if (length(phase_files) == 0) {
    message("PHASE=", phase, ": no DE_workflow files found. Skipping.")
    return(invisible(NULL))
  }

  meta <- tibble(
    file = phase_files,
    dataset_id = map_chr(phase_files, get_dataset_from_path)
  ) %>%
    mutate(is_human = str_detect(dataset_id, "_HUMAN$")) %>%
    filter(
      !is_human,
      dataset_id %in% LOCKED_DATASETS_MOUSE
    )

  if (nrow(meta) == 0) {
    message("PHASE=", phase, ": found files, but none match locked mouse datasets: ",
            paste(LOCKED_DATASETS_MOUSE, collapse = ", "), ". Skipping.")
    return(invisible(NULL))
  }

  # Audit: contrast counts
  by_ds <- meta %>%
    count(dataset_id, name = "n_contrasts") %>%
    arrange(desc(n_contrasts))
  overall <- tibble(
    phase = phase,
    datasets_present = n_distinct(meta$dataset_id),
    total_contrasts = nrow(meta)
  )
  write_tsv(by_ds, file.path(out_root, "contrast_counts_by_dataset.tsv"))
  write_tsv(overall, file.path(out_root, "contrast_counts_overall.tsv"))

  datasets_present <- sort(unique(meta$dataset_id))
  K_present <- length(datasets_present)

  if (K_present < 2) {
    stop("PHASE=", phase, ": fewer than 2 locked datasets present. Need >=2 for consensus.\n",
         "Present: ", paste(datasets_present, collapse = ", "))
  }

  required_support <- ceiling(K_present * consensus_k_fraction)

  message("\nPHASE=", phase,
          " | locked datasets present: ", paste(datasets_present, collapse = ", "),
          " | K_present=", K_present,
          " | core requires >= ", required_support, " supporting datasets")

  # Load DE + pass flag
  de_long <- meta %>%
    mutate(de = map(file, read_de_workflow)) %>%
    tidyr::unnest(de) %>%
    mutate(
      pass = !is.na(FDR) & (FDR <= fdr_cutoff) &
             !is.na(log2FC) & (log2FC >= lfc_cutoff)
    )

  # Within-dataset reproducibility
  support_by_dataset <- de_long %>%
  mutate(pass = ifelse(is.na(pass), FALSE, pass)) %>%   # NA -> FAIL
  group_by(dataset_id, gene_id) %>%
  summarise(
    n_contrasts_dataset = n_distinct(file),
    n_pass = sum(pass),
    support_in_dataset = n_pass / n_contrasts_dataset,  # denominator includes NA contrasts
    dataset_support_flag = support_in_dataset >= within_dataset_support_cutoff,
    .groups = "drop"
  )
  
  # Across-dataset consensus
  core_support <- support_by_dataset %>%
    group_by(gene_id) %>%
    summarise(
      support_datasets = sum(dataset_support_flag, na.rm = TRUE),
      support_fraction = support_datasets / K_present,
      datasets_supporting = paste(sort(unique(dataset_id[dataset_support_flag])), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(support_datasets), desc(support_fraction), gene_id)

  core_gene_set <- core_support %>%
    filter(support_datasets >= required_support) %>%
    arrange(desc(support_datasets), desc(support_fraction), gene_id)

  if (nrow(core_gene_set) == 0) {
    stop("PHASE=", phase, ": core gene set is EMPTY under thresholds.\n",
         "Try lowering lfc_cutoff or within_dataset_support_cutoff.\n",
         "Current: FDR<=", fdr_cutoff,
         ", log2FC>=", lfc_cutoff,
         ", within_dataset_support_cutoff=", within_dataset_support_cutoff,
         ", across>=ceil(K_present*", consensus_k_fraction, ")")
  }

  out_path <- file.path(out_root, "core_gene_set.tsv")
  write_tsv(core_gene_set, out_path)

  message("PHASE=", phase, ": wrote core gene set: ", out_path, " (n=", nrow(core_gene_set), ")")
  invisible(NULL)
}

# -------------------------
# MAIN
# -------------------------
if (!dir.exists(de_comp_root)) {
  stop("DE comparison root not found: ", de_comp_root,
       "\nExpected: <project_root>/04_de/comparison to exist.")
}

for (ph in PHASES) {
  build_core_for_phase(ph)
}

message("\nDONE: Step 06A core gene sets built for phases: ", paste(PHASES, collapse = ", "))
