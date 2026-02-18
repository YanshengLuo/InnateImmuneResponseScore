#!/usr/bin/env Rscript
# ============================================================
# Contrast-level splitting (DESeq2-ready) + optional control fallback
#
# For each input design TSV, generate ONE design file per:
#   tissue + timepoint_hr + batch + delivery_group
# containing ALL delivery replicates + ALL controls.
#
# IMPORTANT:
# Some datasets (e.g., GSE264344) only have baseline/control at timepoint 0.
# If no controls exist at a given timepoint stratum, we FALL BACK to time=0
# controls within the same tissue + batch (configurable).
#
# Windows-safe short paths/filenames to avoid "Filename too long".
# ============================================================

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

# Control fallback mode:
# "strict"  -> require controls in same tissue+time+batch (no fallback)
# "time0"   -> if no controls at that time, use controls at time=0 within same tissue+batch
CONTROL_FALLBACK_MODE <- "time0"

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

# stable short hash (no extra packages)
short_hash <- function(x) {
  ints <- utf8ToInt(as.character(x))
  if (length(ints) == 0) return("00000000")
  h <- sum(ints * seq_along(ints)) %% 100000000
  sprintf("%08d", h)
}

# if group label too long, replace with GID########
short_group_token <- function(g, max_len = 35) {
  g2 <- safe_token(g)
  if (nchar(g2) <= max_len) return(g2)
  paste0("GID", short_hash(g2))
}

# parse numeric hours from timepoint_hr or time_label
to_hours <- function(x) {
  # x can be "72", "72h", "0", "0h", etc.
  x <- as.character(x)
  x <- str_replace(x, "h$", "")
  suppressWarnings(as.numeric(x))
}

# ----------------------------
# Control detection
# ----------------------------
# baseline_0 is your main control when present.
# These two are known controls in your example files.
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

  # Priority 1: baseline_0 if present
  if ("baseline_0" %in% groups) controls <- c(controls, "baseline_0")

  # Priority 2: file-specific known control labels
  if (file_base %in% names(known_controls_by_file)) {
    controls <- c(controls, known_controls_by_file[[file_base]])
  }

  controls <- unique(controls)
  if (length(controls) > 0) return(controls)

  # Priority 3: keyword fallback
  gnorm <- norm_label(groups)
  groups[str_detect(gnorm, control_regex)] %>% unique()
}

# ============================
# MAIN LOOP
# ============================
tsv_files <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__split_index\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__contrast_index\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "__group_map\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "^control_report\\.tsv$")]
tsv_files <- tsv_files[!str_detect(basename(tsv_files), "split_design")]

if (length(tsv_files) == 0) stop("No TSV files found in: ", in_dir)

message("Found ", length(tsv_files), " TSV files in: ", in_dir)
message("Output to: ", out_dir)
message("CONTROL_FALLBACK_MODE = ", CONTROL_FALLBACK_MODE)

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

  # Normalize NA batch to "NA" (critical)
  df2 <- df %>%
    mutate(
      sample_id = as.character(sample_id),
      group_raw  = as.character(group),
      tissue_raw = as.character(tissue),

      time_raw   = as.character(timepoint_hr),
      time_h     = to_hours(time_raw),
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
    mutate(is_control = group_raw %in% control_labels)

  # SHORT dataset output folder
  dataset_folder <- file.path(out_dir, safe_token(stem))
  dir.create(dataset_folder, showWarnings = FALSE, recursive = TRUE)

  written <- character(0)
  group_map_rows <- list()

  # Strata: tissue + time + batch (we write one contrast file per delivery group within each stratum)
  strata <- df2 %>%
    distinct(tissue_raw, time_label, batch_raw) %>%
    arrange(tissue_raw, time_label, batch_raw)

  for (s in seq_len(nrow(strata))) {
    tissue_s <- strata$tissue_raw[s]
    time_s   <- strata$time_label[s]
    batch_s  <- strata$batch_raw[s]

    block <- df2 %>%
      filter(tissue_raw == tissue_s, time_label == time_s, batch_raw == batch_s)

    # Controls matched to this stratum
    ctrl_block <- block %>% filter(is_control)

    # FALLBACK: if no controls at this timepoint, use time=0 controls within same tissue+batch
    control_time_policy <- "matched_time"
    if (nrow(ctrl_block) == 0 && CONTROL_FALLBACK_MODE == "time0") {
      ctrl_block <- df2 %>%
        filter(
          tissue_raw == tissue_s,
          batch_raw  == batch_s,
          is_control,
          !is.na(time_h),
          time_h == 0
        )
      if (nrow(ctrl_block) > 0) control_time_policy <- "fallback_time0"
    }

    if (nrow(ctrl_block) == 0) next

    # Delivery groups in this stratum
    delivery_groups <- block %>%
      filter(!is_control) %>%
      pull(group_raw) %>%
      unique() %>%
      sort()

    if (length(delivery_groups) == 0) next

    # Create ONE file per delivery group (full contrast)
    for (g in delivery_groups) {
      deliv_block <- block %>% filter(!is_control, group_raw == g)
      if (nrow(deliv_block) == 0) next

      # Full contrast design: all controls + all delivery reps
      sub <- bind_rows(ctrl_block, deliv_block) %>%
        mutate(
          condition_simple = ifelse(is_control, "CONTROL", "DELIVERY"),
          control_label    = control_labels[1],
          contrast_label   = paste0(g, "__vs__", control_labels[1]),
          control_time_policy = control_time_policy
        )

      # Tokens for naming
      tissue_tok <- safe_token(tissue_s)
      # store hours as number if possible to shorten
      time_tok_num <- to_hours(time_s)
      time_tok <- ifelse(!is.na(time_tok_num), as.character(time_tok_num), safe_token(time_s))
      batch_tok <- safe_token(batch_s)

      g_tok <- short_group_token(g)
      ctrl_tok <- safe_token(control_labels[1])

      # Group mapping
      group_map_rows[[length(group_map_rows) + 1]] <- tibble(
        tissue_raw = tissue_s,
        time_label = time_s,
        batch_raw  = batch_s,
        group_raw  = g,
        group_token = g_tok
      )

      # SHORT filename
      out_name <- sprintf(
        "%s__T=%s__H=%s__B=%s__G=%s__VS=%s.tsv",
        safe_token(stem), tissue_tok, time_tok, batch_tok, g_tok, ctrl_tok
      )
      out_path <- file.path(dataset_folder, out_name)

      write_tsv(sub, out_path)
      written <- c(written, out_path)
    }
  }

  if (length(written) > 0) {
    idx_path <- file.path(dataset_folder, paste0(safe_token(stem), "__contrast_index.tsv"))
    map_path <- file.path(dataset_folder, paste0(safe_token(stem), "__group_map.tsv"))

    gm <- bind_rows(group_map_rows) %>%
      distinct() %>%
      arrange(tissue_raw, time_label, batch_raw, group_token, group_raw)

    idx <- tibble(design_file = written) %>%
      mutate(
        file_base = basename(design_file),
        tissue = str_match(file_base, "__T=([^_]+)__H=")[,2],
        time_h = str_match(file_base, "__H=([^_]+)__B=")[,2],
        batch = str_match(file_base, "__B=([^_]+)__G=")[,2],
        group_token = str_match(file_base, "__G=([^_]+)__VS=")[,2],
        control_token = str_match(file_base, "__VS=([^\\.]+)\\.tsv$")[,2]
      ) %>%
      left_join(gm %>% select(group_token, group_raw),
                by = c("group_token" = "group_token")) %>%
      mutate(group_raw = ifelse(is.na(group_raw), group_token, group_raw))

    write_tsv(idx, idx_path)
    write_tsv(gm, map_path)

    message("Wrote ", length(written), " contrast-level design files.")
    message("Index: ", idx_path)
    message("Group map: ", map_path)
  } else {
    message("No contrast files written for: ", file_base,
            " (controls exist but no delivery groups within strata or no eligible fallback controls).")
  }
}

# Control report across all datasets
report_path <- file.path(out_dir, "control_report.tsv")
tibble(
  file = names(control_report),
  controls = vapply(control_report, function(x) paste(x$control_labels, collapse = ";"), character(1))
) %>% write_tsv(report_path)

message("\nAll done.")
message("Control report: ", report_path)
