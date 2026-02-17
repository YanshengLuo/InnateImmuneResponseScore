#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================
# USER SETTINGS (Windows)
# ============================
in_dir  <- "C:/Users/john/Desktop/IMRS_Project/Hypergator_scripts/InnateImmuneResponseScore/verified_metadata"
out_dir <- file.path(in_dir, "splited")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
norm_label <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

safe_token <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- str_trim(x)
  x <- str_replace_all(x, "[\\\\/:*?\"<>|]+", "_")
  x <- str_replace_all(x, "\\s+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_+|_+$", "")
  ifelse(nchar(x) == 0, "NA", x)
}

pair_controls <- function(deliv_df, ctrl_df) {
  if (nrow(deliv_df) == 0 || nrow(ctrl_df) == 0) return(list())
  deliv_df <- deliv_df %>% arrange(sample_id)
  ctrl_df  <- ctrl_df  %>% arrange(sample_id)

  pairs <- vector("list", nrow(deliv_df))
  for (i in seq_len(nrow(deliv_df))) {
    d <- deliv_df[i, , drop = FALSE]
    cidx <- ((i - 1) %% nrow(ctrl_df)) + 1
    c <- ctrl_df[cidx, , drop = FALSE]
    pairs[[i]] <- list(delivery = d, control = c, k = i)
  }
  pairs
}

# ----------------------------
# Control detection rules
# ----------------------------
# 1) baseline_0 always control if present
# 2) otherwise, per-file known controls (from your datasets)
# 3) otherwise keyword/regex fallback
known_controls_by_file <- list(
  "GSE190850_HUMAN_design.tsv" = c("delivery_unstimulated"),
  "GSE39129_design.tsv"        = c("delivery_none")
)

control_keywords <- c(
  "baseline", "control", "ctrl",
  "pbs", "saline", "vehicle",
  "mock", "sham", "empty",
  "naive", "untreated", "uninfected",
  "unstimulated", "none", "no_treatment", "pre"
)
control_regex <- paste0("(", paste(unique(norm_label(control_keywords)), collapse = "|"), ")")

detect_controls <- function(df, file_base) {
  groups <- unique(as.character(df$group))
  controls <- character(0)

  # priority 1
  if ("baseline_0" %in% groups) controls <- c(controls, "baseline_0")

  # priority 2: known per file
  if (file_base %in% names(known_controls_by_file)) {
    controls <- c(controls, known_controls_by_file[[file_base]])
  }

  controls <- unique(controls)
  if (length(controls) > 0) return(controls)

  # priority 3: regex fallback
  gnorm <- norm_label(groups)
  controls <- groups[str_detect(gnorm, control_regex)]
  unique(controls)
}

# ============================
# MAIN LOOP
# ============================
tsv_files <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__split_index\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "split_design")]

if (length(tsv_files) == 0) stop("No TSV files found in: ", in_dir)

message("Found ", length(tsv_files), " TSV files in: ", in_dir)
message("Output to: ", out_dir)

# optional: write a quick control report
control_report <- list()

for (f in tsv_files) {
  file_base <- basename(f)
  stem <- tools::file_path_sans_ext(file_base)
  message("\n=== Processing: ", file_base, " ===")

  df <- read_tsv(f, show_col_types = FALSE)
  req <- c("sample_id", "group", "tissue", "timepoint_hr", "batch")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    message("SKIP (missing columns): ", paste(miss, collapse = ", "))
    next
  }

  # FIX: normalize NA batch to string so stratification/filtering works
  df2 <- df %>%
    mutate(
      group_raw  = as.character(group),
      tissue_raw = as.character(tissue),
      time_raw   = as.character(timepoint_hr),
      time_h     = suppressWarnings(as.numeric(time_raw)),
      time_label = ifelse(!is.na(time_h), paste0(time_h, "h"), safe_token(time_raw)),
      batch_raw  = as.character(batch),
      batch_raw  = ifelse(is.na(batch_raw) | batch_raw == "", "NA", batch_raw)
    )

  control_labels <- detect_controls(df2, file_base)
  control_report[[file_base]] <- list(control_labels = control_labels)

  if (length(control_labels) == 0) {
    message("WARN: No controls detected in ", file_base, ". Skipping.")
    next
  } else {
    message("Controls detected: ", paste(control_labels, collapse = ", "))
  }

  df2 <- df2 %>%
    mutate(
      is_control = group_raw %in% control_labels
    )

  strata <- df2 %>%
    distinct(tissue_raw, time_label, batch_raw) %>%
    arrange(tissue_raw, time_label, batch_raw)

  written <- character(0)

  for (s in seq_len(nrow(strata))) {
    tissue_s <- strata$tissue_raw[s]
    time_s   <- strata$time_label[s]
    batch_s  <- strata$batch_raw[s]

    block <- df2 %>%
      filter(tissue_raw == tissue_s, time_label == time_s, batch_raw == batch_s)

    ctrl_block <- block %>% filter(is_control)
    if (nrow(ctrl_block) == 0) next

    delivery_groups <- block %>%
      filter(!is_control) %>%
      pull(group_raw) %>%
      unique() %>%
      sort()

    if (length(delivery_groups) == 0) next

    for (g in delivery_groups) {
      deliv_block <- block %>% filter(group_raw == g, !is_control)
      pairs <- pair_controls(deliv_block, ctrl_block)
      if (length(pairs) == 0) next

      tissue_tok <- safe_token(tissue_s)
      time_tok   <- safe_token(time_s)
      batch_tok  <- safe_token(batch_s)
      g_tok      <- safe_token(g)

      # For naming, use the FIRST control label as "primary"
      ctrl_name_primary <- safe_token(control_labels[1])

      contrast_folder <- file.path(
        out_dir,
        safe_token(stem),
        paste0("tissue=", tissue_tok),
        paste0("time=", time_tok),
        paste0("batch=", batch_tok),
        paste0("contrast=", g_tok, "__vs__", ctrl_name_primary)
      )
      dir.create(contrast_folder, showWarnings = FALSE, recursive = TRUE)

      for (p in pairs) {
        d <- p$delivery
        c <- p$control
        k <- p$k

        sub <- bind_rows(c, d) %>%
          mutate(
            condition_simple = ifelse(is_control, "CONTROL", "DELIVERY"),
            control_label    = control_labels[1],
            contrast_label   = paste0(g, "__vs__", control_labels[1]),
            paired_delivery_sample = d$sample_id,
            paired_control_sample  = c$sample_id
          )

        out_name <- paste0(
          safe_token(stem),
          "__tissue=", tissue_tok,
          "__time=", time_tok,
          "__batch=", batch_tok,
          "__contrast=", g_tok, "__vs__", ctrl_name_primary,
          "__pair=", k,
          "__D=", safe_token(d$sample_id),
          "__C=", safe_token(c$sample_id),
          ".tsv"
        )

        out_path <- file.path(contrast_folder, out_name)
        write_tsv(sub, out_path)
        written <- c(written, out_path)
      }
    }
  }

  if (length(written) > 0) {
    idx_path <- file.path(out_dir, safe_token(stem), paste0(safe_token(stem), "__split_index.tsv"))
    dir.create(dirname(idx_path), showWarnings = FALSE, recursive = TRUE)
    tibble(split_file = written) %>% write_tsv(idx_path)
    message("Wrote ", length(written), " split files. Index: ", idx_path)
  } else {
    message("No split files written for: ", file_base, " (controls exist but maybe no delivery groups within strata).")
  }
}

# write a simple control report
report_path <- file.path(out_dir, "control_report.tsv")
tibble(
  file = names(control_report),
  controls = vapply(control_report, function(x) paste(x$control_labels, collapse = ";"), character(1))
) %>% write_tsv(report_path)

message("\nAll done.")
message("Control report: ", report_path)
