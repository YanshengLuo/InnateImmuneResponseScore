suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# =========================
# Config
# =========================
imrs_root <- "D:/IMRS_Project"
framework_root <- "D:/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/Dataset_Level_Confidence_score_prototype"

de_root <- file.path(imrs_root, "04_de", "comparison")
manifest_file <- file.path(framework_root, "inputs", "01_dataset_manifest.tsv")

out_file <- file.path(framework_root, "inputs", "05_gene_effects.tsv")
log_file <- file.path(framework_root, "logs", "04_build_gene_effects.log")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

# =========================
# Logging helper
# =========================
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

# =========================
# Helpers
# =========================
safe_read_table <- function(path) {
  ext <- tools::file_ext(path)
  if (tolower(ext) == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  } else {
    return(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  }
}

pick_first_col <- function(df, candidates) {
  hits <- names(df)[tolower(names(df)) %in% tolower(candidates)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

score_candidate_file <- function(path) {
  b <- tolower(basename(path))
  score <- 0
  if (grepl("workflow", b)) score <- score + 5
  if (grepl("schema", b)) score <- score + 5
  if (grepl("delivery", b)) score <- score + 2
  if (grepl("control", b)) score <- score + 2
  if (grepl("contrast", b)) score <- score + 1
  if (grepl("full", b)) score <- score - 1
  score
}

# =========================
# Inputs
# =========================
if (!file.exists(manifest_file)) {
  stop("Missing manifest file: ", manifest_file)
}
if (!dir.exists(de_root)) {
  stop("Missing DE root: ", de_root)
}

manifest <- read_tsv(manifest_file, show_col_types = FALSE, progress = FALSE)

dataset_ids <- manifest %>%
  filter(isTRUE(include_in_analysis) | include_in_analysis == TRUE | include_in_analysis == "TRUE") %>%
  pull(dataset_id) %>%
  unique()

log_message("DE root: ", de_root)
log_message("Datasets to scan: ", length(dataset_ids))

all_rows <- list()

# =========================
# Main loop
# =========================
for (dataset_id in dataset_ids) {
  ds_dir <- file.path(de_root, dataset_id)

  if (!dir.exists(ds_dir)) {
    log_message("Missing DE directory for ", dataset_id, ": ", ds_dir)
    next
  }

  candidate_files <- list.files(
    ds_dir,
    pattern = "\\.(tsv|txt|csv)$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(candidate_files) == 0) {
    log_message("No candidate DE files found for ", dataset_id)
    next
  }

  candidate_tbl <- tibble(
    file = candidate_files,
    score = vapply(candidate_files, score_candidate_file, numeric(1))
  ) %>%
    arrange(desc(score), file)

  selected_df <- NULL
  selected_file <- NA_character_

  for (f in candidate_tbl$file) {
    df <- tryCatch(safe_read_table(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) next

    gene_col <- pick_first_col(df, c("gene", "gene_id"))
    effect_col <- pick_first_col(df, c("log2FC", "log2fc", "beta", "effect"))
    padj_col <- pick_first_col(df, c("padj", "FDR", "fdr"))
    stat_col <- pick_first_col(df, c("stat", "z", "t", "wald_stat"))

    if (!is.na(gene_col) && !is.na(effect_col)) {
      selected_df <- df
      selected_file <- f
      log_message("Selected DE file for ", dataset_id, ": ", f)
      break
    }
  }

  if (is.null(selected_df)) {
    log_message("No usable DE table found for ", dataset_id)
    next
  }

  gene_col <- pick_first_col(selected_df, c("gene", "gene_id"))
  effect_col <- pick_first_col(selected_df, c("log2FC", "log2fc", "beta", "effect"))
  padj_col <- pick_first_col(selected_df, c("padj", "FDR", "fdr"))
  stat_col <- pick_first_col(selected_df, c("stat", "z", "t", "wald_stat"))

  out_df <- selected_df %>%
    transmute(
      dataset_id = dataset_id,
      gene_id = as.character(.data[[gene_col]]),
      effect = suppressWarnings(as.numeric(.data[[effect_col]])),
      effect_type = "log2FC",
      padj = if (!is.na(padj_col)) suppressWarnings(as.numeric(.data[[padj_col]])) else NA_real_,
      stat = if (!is.na(stat_col)) suppressWarnings(as.numeric(.data[[stat_col]])) else NA_real_,
      source_file = selected_file
    ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    filter(!is.na(effect))

  log_message("Rows kept for ", dataset_id, ": ", nrow(out_df))
  all_rows[[length(all_rows) + 1]] <- out_df
}

gene_effects <- bind_rows(all_rows)

if (nrow(gene_effects) == 0) {
  stop("No gene effects were successfully built.")
}

write_tsv(gene_effects, out_file)
log_message("Saved gene effects: ", out_file)
log_message("Rows written: ", nrow(gene_effects))
log_message("Datasets represented: ", n_distinct(gene_effects$dataset_id))