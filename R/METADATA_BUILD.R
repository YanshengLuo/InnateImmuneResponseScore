## ================================
## METADATA_BUILD.R
## 00_metadata construction (robust)
##
## Dataset: <dataset_id>
##
## Reads:
##   00_metadata/<dataset_id>/SraRunTable.csv
##
## Writes:
##   00_metadata/<dataset_id>/<dataset_id>_samples.tsv
##   00_metadata/<dataset_id>/<dataset_id>_design.tsv
##
## Key logic:
##   1. Detect SRR/BioSample/Experiment columns.
##   2. If treatment/delivery column exists, use it.
##   3. Also build metadata_blob from ALL columns.
##   4. Detect CONTROL from dictionary terms anywhere in metadata_blob.
##   5. Otherwise assign DELIVERY.
##   6. Write condition_simple for downstream DE/scoring.
## ================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

## ----------------------------
## 1) Dataset + paths
## ----------------------------
if (!exists("dataset_id", inherits = FALSE)) {
  stop("dataset_id is not set. Run via wrapper or set dataset_id before sourcing METADATA_BUILD.R.")
}

project_root <- "D:/IMRS_Project"
meta_dir <- file.path(project_root, "00_metadata", dataset_id)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

runinfo_path <- file.path(meta_dir, "SraRunTable.csv")
if (!file.exists(runinfo_path)) {
  stop("RunInfo file not found: ", runinfo_path)
}

## ----------------------------
## 2) Helpers
## ----------------------------
pick_col <- function(df, candidates) {
  nms <- names(df)
  dn  <- tolower(nms)
  cc  <- tolower(candidates)

  hit <- intersect(cc, dn)

  if (length(hit) == 0) {
    return(NA_character_)
  }

  nms[match(hit[1], dn)]
}

report_missing <- function(col_map, required, context = "") {
  cat("\n================ METADATA COLUMN DETECTION REPORT ================\n")

  if (nzchar(context)) {
    cat("Context:", context, "\n\n")
  }

  for (nm in names(col_map)) {
    cat(sprintf(
      "  %-14s -> %s\n",
      nm,
      ifelse(is.na(col_map[[nm]]), "NOT FOUND", col_map[[nm]])
    ))
  }

  missing_required <- required[is.na(unlist(col_map[required]))]

  cat("\nRequired fields:\n")

  if (length(missing_required) == 0) {
    cat("  OK: all required fields detected.\n")
  } else {
    cat("  MISSING required:", paste(missing_required, collapse = ", "), "\n")
  }

  cat("===============================================================\n\n")

  missing_required
}

assert_nonempty <- function(df, col) {
  bad <- is.na(df[[col]]) | trimws(as.character(df[[col]])) == ""

  if (any(bad)) {
    stop(sprintf(
      "Column '%s' has %d missing/empty values.",
      col,
      sum(bad)
    ))
  }
}

make_metadata_blob <- function(df) {
  df_chr <- df %>%
    mutate(across(everything(), as.character))

  apply(df_chr, 1, function(x) {
    x <- x[!is.na(x)]
    x <- paste(x, collapse = " | ")
    x <- str_to_lower(x)
    x <- str_replace_all(x, "\u00A0", " ")
    x <- str_replace_all(x, "[^a-z0-9\\.]+", " ")
    str_squish(x)
  })
}

parse_timepoint_hr <- function(timepoint_raw) {
  tp0 <- ifelse(is.na(timepoint_raw), "", as.character(timepoint_raw))
  tp0 <- str_trim(str_replace_all(tp0, "\u00A0", " "))

  tp <- str_to_lower(tp0)
  tp <- str_replace_all(tp, "[^a-z0-9\\.]+", " ")
  tp <- str_squish(tp)

  is_zero <- tp %in% c(
    "baseline",
    "base line",
    "pre",
    "pretreatment",
    "pre treatment",
    "pre vaccination",
    "prevaccination",
    "week 0",
    "wk 0",
    "day 0",
    "naive",
    "untreated",
    "control",
    "vehicle",
    "pbs",
    "saline",
    "nothing",
    "nothing injected",
    "mock",
    "sham",
    "uninfected"
  )

  out <- rep(NA_real_, length(tp))
  out[is_zero] <- 0

  need <- !is_zero & tp != ""

  has_hr_unit <- need & str_detect(
    tp,
    "^[0-9]+\\.?[0-9]*\\s*(h|hr|hrs|hour|hours)\\b|\\b[0-9]+\\.?[0-9]*\\s*(h|hr|hrs|hour|hours)\\b"
  )

  out[has_hr_unit] <- suppressWarnings(parse_number(tp[has_hr_unit]))

  is_numeric_only <- need & str_detect(tp, "^[0-9]+\\.?[0-9]*$")

  out[is_numeric_only] <- suppressWarnings(parse_number(tp[is_numeric_only]))

  out
}

detect_control_from_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- str_to_lower(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[^a-z0-9\\.]+", " ")
  x <- str_squish(x)

  str_detect(x, "\\bbaseline\\b") |
    str_detect(x, "\\bpre\\b|\\bpre vaccination\\b|\\bprevaccination\\b") |
    str_detect(x, "\\bday 0\\b|\\bweek 0\\b|\\bwk 0\\b") |
    str_detect(x, "\\buntreated\\b|\\bnaive\\b|\\bnormal\\b") |
    str_detect(x, "\\bpbs\\b|\\bsaline\\b|\\bvehicle\\b|\\bmock\\b|\\bsham\\b") |
    str_detect(x, "\\bnothing\\b|\\bnothing injected\\b|\\buninfected\\b") |
    str_detect(x, "\\bcontrol\\b")
}

detect_time_from_blob <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- str_to_lower(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[^a-z0-9\\.]+", " ")
  x <- str_squish(x)

  out <- rep(NA_real_, length(x))

  zero_flag <- detect_control_from_text(x)
  out[zero_flag] <- 0

  needs <- is.na(out)

  hour_hit <- str_extract(
    x,
    "\\b[0-9]+\\.?[0-9]*\\s*(h|hr|hrs|hour|hours)\\b"
  )

  out[needs & !is.na(hour_hit)] <- suppressWarnings(parse_number(hour_hit[needs & !is.na(hour_hit)]))

  out
}

## ----------------------------
## 3) Read RunInfo
## ----------------------------
runinfo <- read_csv(runinfo_path, show_col_types = FALSE)

cat("\nLoaded RunInfo:\n")
cat("  File:", runinfo_path, "\n")
cat("  Rows:", nrow(runinfo), "  Cols:", ncol(runinfo), "\n")

if (nrow(runinfo) == 0) {
  stop("SraRunTable.csv has zero rows: ", runinfo_path)
}

runinfo$metadata_blob <- make_metadata_blob(runinfo)

## ----------------------------
## 4) Column synonyms
## ----------------------------
syn <- list(
  srr_id = c(
    "Run",
    "run",
    "SRR",
    "srr",
    "srr_id",
    "srr_accession"
  ),

  biosample_id = c(
    "BioSample",
    "biosample",
    "BioSample_ID",
    "biosample_id"
  ),

  experiment_id = c(
    "Experiment",
    "experiment",
    "SRX",
    "srx",
    "experiment_id"
  ),

  tissue = c(
    "tissue",
    "organ",
    "Organ",
    "tissue_type",
    "source_name",
    "source_name_ch1"
  ),

  source_name = c(
    "source_name",
    "source_name_ch1",
    "SourceName",
    "source",
    "Sample Name",
    "sample_name",
    "title"
  ),

  timepoint = c(
    "timepoint",
    "time_point",
    "time",
    "Time",
    "collection_time",
    "hours_post",
    "hrs_post",
    "time after treatment",
    "time_after_treatment",
    "treatment_time"
  ),

  delivery_type = c(
    "vaccine_type",
    "vaccine",
    "treatment",
    "treat",
    "platform",
    "vector",
    "delivery",
    "group",
    "arm",
    "condition",
    "sample group",
    "sample_group",
    "title",
    "Sample Name",
    "sample_name",
    "source_name",
    "source_name_ch1"
  )
)

col_map <- list(
  srr_id        = pick_col(runinfo, syn$srr_id),
  biosample_id  = pick_col(runinfo, syn$biosample_id),
  experiment_id = pick_col(runinfo, syn$experiment_id),
  tissue        = pick_col(runinfo, syn$tissue),
  source_name   = pick_col(runinfo, syn$source_name),
  timepoint     = pick_col(runinfo, syn$timepoint),
  delivery_type = pick_col(runinfo, syn$delivery_type)
)

required_detect <- c("srr_id", "biosample_id", "experiment_id")

missing_required <- report_missing(
  col_map,
  required_detect,
  context = paste(dataset_id, "RunInfo mapping")
)

if (length(missing_required) > 0) {
  stop(
    "Stopping: required metadata fields were not detected: ",
    paste(missing_required, collapse = ", ")
  )
}

timepoint_missing <- is.na(col_map$timepoint)
delivery_type_missing <- is.na(col_map$delivery_type)

## ----------------------------
## 5) Build standardized samples table
## ----------------------------
samples <- runinfo %>%
  transmute(
    dataset_id = dataset_id,

    srr_id = as.character(.data[[col_map$srr_id]]),
    biosample_id = as.character(.data[[col_map$biosample_id]]),
    experiment_id = as.character(.data[[col_map$experiment_id]]),

    metadata_blob = metadata_blob,

    species = "Mus_musculus",

    tissue = if (!is.na(col_map$tissue)) {
      as.character(.data[[col_map$tissue]])
    } else {
      NA_character_
    },

    source_name = if (!is.na(col_map$source_name)) {
      as.character(.data[[col_map$source_name]])
    } else {
      NA_character_
    },

    timepoint_raw = if (!timepoint_missing) {
      as.character(.data[[col_map$timepoint]])
    } else {
      NA_character_
    },

    delivery_type = if (!delivery_type_missing) {
      as.character(.data[[col_map$delivery_type]])
    } else {
      NA_character_
    },

    condition = NA_character_,
    control_type = NA_character_,
    payload = NA_character_,

    replicate = NA_integer_,
    batch = NA_character_,
    notes = NA_character_
  )

assert_nonempty(samples, "srr_id")

if (!delivery_type_missing) {
  # Do not stop on empty delivery_type because metadata_blob can rescue it.
  empty_delivery <- is.na(samples$delivery_type) | trimws(samples$delivery_type) == ""
  if (any(empty_delivery)) {
    warning(
      "delivery_type column detected, but ",
      sum(empty_delivery),
      " rows are empty. metadata_blob rescue will be used."
    )
  }
}

## ----------------------------
## 5B) Optional manual annotation override
## ----------------------------
manual_path <- file.path(meta_dir, paste0(dataset_id, "_manual_annotations.tsv"))

if (file.exists(manual_path)) {
  cat("\nManual annotation file detected:\n")
  cat("  ", manual_path, "\n")

  manual <- read_tsv(manual_path, show_col_types = FALSE)

  required_manual <- c("srr_id", "condition_simple")
  missing_manual <- setdiff(required_manual, names(manual))

  if (length(missing_manual) > 0) {
    stop(
      "Manual annotation file is missing required columns: ",
      paste(missing_manual, collapse = ", ")
    )
  }

  manual <- manual %>%
    mutate(
      srr_id = as.character(srr_id),
      condition_simple = toupper(trimws(as.character(condition_simple))),

      delivery_type = if ("delivery_type" %in% names(.)) {
        as.character(delivery_type)
      } else {
        NA_character_
      },

      tissue = if ("tissue" %in% names(.)) {
        as.character(tissue)
      } else {
        NA_character_
      },

      timepoint_hr_manual = if ("timepoint_hr" %in% names(.)) {
        suppressWarnings(as.numeric(timepoint_hr))
      } else {
        NA_real_
      },

      batch_manual = if ("batch" %in% names(.)) {
        as.character(batch)
      } else {
        NA_character_
      },

      notes_manual = if ("notes" %in% names(.)) {
        as.character(notes)
      } else {
        NA_character_
      }
    )

  bad_condition <- manual %>%
    filter(!condition_simple %in% c("CONTROL", "DELIVERY"))

  if (nrow(bad_condition) > 0) {
    stop("Manual annotation condition_simple must be CONTROL or DELIVERY only.")
  }

  samples <- samples %>%
    left_join(
      manual %>%
        select(
          srr_id,
          condition_simple_manual = condition_simple,
          delivery_type_manual = delivery_type,
          tissue_manual = tissue,
          timepoint_hr_manual,
          batch_manual,
          notes_manual
        ),
      by = "srr_id"
    ) %>%
    mutate(
      condition = case_when(
        !is.na(condition_simple_manual) & condition_simple_manual == "CONTROL" ~ "control",
        !is.na(condition_simple_manual) & condition_simple_manual == "DELIVERY" ~ "delivery",
        TRUE ~ condition
      ),

      delivery_type = coalesce(delivery_type_manual, delivery_type),
      tissue = coalesce(tissue_manual, tissue),
      batch = coalesce(batch_manual, batch),

      notes = case_when(
        !is.na(notes_manual) & notes_manual != "" ~ notes_manual,
        TRUE ~ notes
      )
    ) %>%
    select(
      -condition_simple_manual,
      -delivery_type_manual,
      -tissue_manual,
      -batch_manual,
      -notes_manual
    )
} else {
  samples$timepoint_hr_manual <- NA_real_
}

## ----------------------------
## 6) Parse timepoints → hours
## ----------------------------
samples <- samples %>%
  mutate(
    delivery_type = ifelse(is.na(delivery_type), "", delivery_type),
    delivery_type = str_squish(str_trim(delivery_type)),

    metadata_blob = ifelse(is.na(metadata_blob), "", metadata_blob),
    metadata_blob = str_to_lower(metadata_blob),

    timepoint_raw = ifelse(
      is.na(timepoint_raw),
      NA_character_,
      str_trim(str_replace_all(timepoint_raw, "\u00A0", " "))
    ),

    timepoint_hr_from_column = if (timepoint_missing) {
      NA_real_
    } else {
      parse_timepoint_hr(timepoint_raw)
    },

    timepoint_hr_from_blob = detect_time_from_blob(metadata_blob),

    timepoint_hr = coalesce(
      timepoint_hr_manual,
      timepoint_hr_from_column,
      timepoint_hr_from_blob
    ),

    notes = case_when(
      !is.na(timepoint_hr_manual) ~ notes,
      timepoint_missing & is.na(timepoint_hr_from_blob) ~ if_else(
        is.na(notes) | notes == "",
        "timepoint_missing_not_resolved",
        paste(notes, "timepoint_missing_not_resolved", sep = ";")
      ),
      timepoint_missing & !is.na(timepoint_hr_from_blob) ~ if_else(
        is.na(notes) | notes == "",
        "timepoint_missing_rescued_from_metadata_blob",
        paste(notes, "timepoint_missing_rescued_from_metadata_blob", sep = ";")
      ),
      TRUE ~ notes
    )
  )

cat("\nTimepoint diagnostics raw unique:\n")
print(sort(unique(samples$timepoint_raw)))

na_after_parse <- samples %>%
  filter(is.na(timepoint_hr)) %>%
  count(timepoint_raw, sort = TRUE)

if (nrow(na_after_parse) > 0) {
  cat("\nWARNING: timepoint_hr is NA after parsing:\n")
  print(na_after_parse)

  samples <- samples %>%
    mutate(
      notes = if_else(
        is.na(timepoint_hr),
        if_else(
          is.na(notes) | notes == "",
          "timepoint_hr_unresolved",
          paste(notes, "timepoint_hr_unresolved", sep = ";")
        ),
        notes
      )
    )
}

## ----------------------------
## 7) Define condition and group
## ----------------------------
samples <- samples %>%
  mutate(
    delivery_lc = str_to_lower(delivery_type),

    is_control_delivery_col = detect_control_from_text(delivery_lc),
    is_control_blob = detect_control_from_text(metadata_blob),

    condition = case_when(
      !is.na(condition) & condition %in% c("control", "delivery") ~ condition,

      is_control_delivery_col ~ "control",
      is_control_blob ~ "control",

      TRUE ~ "delivery"
    ),

    control_type = if_else(
      condition == "control",
      "baseline_or_vehicle",
      NA_character_
    ),

    group = case_when(
      condition == "control" ~ "baseline_0",

      condition == "delivery" & delivery_lc != "" ~ paste0(
        "delivery_",
        str_replace_all(delivery_lc, "[^a-z0-9]+", "_")
      ),

      condition == "delivery" & delivery_lc == "" ~ "delivery_unknown",

      TRUE ~ "needs_review"
    ),

    notes = case_when(
      condition == "control" & is_control_blob & !is_control_delivery_col ~ if_else(
        is.na(notes) | notes == "",
        "control_rescued_from_metadata_blob",
        paste(notes, "control_rescued_from_metadata_blob", sep = ";")
      ),

      condition == "delivery" & delivery_lc == "" ~ if_else(
        is.na(notes) | notes == "",
        "delivery_type_missing_assigned_delivery",
        paste(notes, "delivery_type_missing_assigned_delivery", sep = ";")
      ),

      TRUE ~ notes
    )
  ) %>%
  select(
    -delivery_lc,
    -is_control_delivery_col,
    -is_control_blob,
    -timepoint_hr_manual,
    -timepoint_hr_from_column,
    -timepoint_hr_from_blob
  )

cat("\nCondition assignment counts:\n")
print(table(samples$condition, useNA = "ifany"))

cat("\nGroup counts:\n")
print(table(samples$group, useNA = "ifany"))

## ----------------------------
## 8) Final checks
## ----------------------------
cat("\nBaseline controls by tissue/timepoint:\n")
print(
  samples %>%
    filter(condition == "control") %>%
    count(tissue, timepoint_hr, sort = TRUE)
)

cat("\nDelivery samples by tissue/delivery/timepoint:\n")
print(
  samples %>%
    filter(condition == "delivery") %>%
    count(tissue, delivery_type, timepoint_hr, sort = TRUE)
)

if (!any(samples$condition == "control")) {
  warning("No control samples detected. DE contrasts may fail for this dataset.")

  samples <- samples %>%
    mutate(
      notes = if_else(
        is.na(notes) | notes == "",
        "no_controls_detected",
        paste(notes, "no_controls_detected", sep = ";")
      )
    )
}

if (!any(samples$condition == "delivery")) {
  warning("No delivery samples detected. DE contrasts may fail for this dataset.")

  samples <- samples %>%
    mutate(
      notes = if_else(
        is.na(notes) | notes == "",
        "no_delivery_detected",
        paste(notes, "no_delivery_detected", sep = ";")
      )
    )
}

## ----------------------------
## 9) Save outputs
## ----------------------------
samples_path <- file.path(meta_dir, paste0(dataset_id, "_samples.tsv"))

samples_out <- samples %>%
  mutate(
    condition_simple = case_when(
      condition == "control" ~ "CONTROL",
      condition == "delivery" ~ "DELIVERY",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    dataset_id,
    srr_id,
    biosample_id,
    experiment_id,
    species,
    tissue,
    source_name,
    timepoint_raw,
    timepoint_hr,
    delivery_type,
    condition,
    condition_simple,
    control_type,
    payload,
    group,
    replicate,
    batch,
    notes,
    metadata_blob
  )

write_tsv(samples_out, samples_path)

design_out <- samples_out %>%
  transmute(
    sample_id = srr_id,
    condition_simple = condition_simple,
    group = group,
    tissue = tissue,
    timepoint_hr = timepoint_hr,
    batch = batch
  )

design_path <- file.path(meta_dir, paste0(dataset_id, "_design.tsv"))
write_tsv(design_out, design_path)

cat("\nWrote:\n")
cat("  -", samples_path, "\n")
cat("  -", design_path, "\n\n")

cat("Done building metadata for:", dataset_id, "\n")