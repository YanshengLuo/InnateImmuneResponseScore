## ================================
## 00_metadata construction
## Dataset: GSE264344
## ================================

library(readr)
library(dplyr)

## ----------------
## 1. Read RunInfo
## ----------------
runinfo <- read_csv(
  "C:/Users/john/Desktop/IMRS_Project/R_codes/data/SraRunTable.csv",
  show_col_types = FALSE
)

# quick sanity checks
colnames(runinfo)
nrow(runinfo)
head(runinfo, 3)

## ----------------
## 2. Initialize dataset
## ----------------
dataset_id <- "GSE264344"

samples <- runinfo %>%
  transmute(
    dataset_id = dataset_id,
    srr_id = Run,
    biosample_id = BioSample,
    experiment_id = Experiment,
    
    species = "Mus_musculus",   # confirmed from paper
    
    tissue = tissue,
    source_name = source_name,
    
    timepoint_raw = timepoint,      # keep raw first
    delivery_type = vaccine_type,   # Ad26 / Ad5 / baseline
    
    condition = NA_character_,      # filled later
    control_type = NA_character_,   # filled later
    payload = NA_character_,        # not needed for this dataset
    
    replicate = NA_integer_,
    batch = NA_character_,
    notes = NA_character_
  )

## ----------------
## 3. Parse timepoints (→ hours)
## ----------------
samples <- samples %>%
  mutate(
    timepoint_raw = trimws(timepoint_raw),
    timepoint_hr = case_when(
      tolower(timepoint_raw) == "baseline" ~ 0,
      grepl("\\bhr\\b|\\bhrs\\b|\\bhour\\b|\\bhours\\b", tolower(timepoint_raw)) ~
        as.numeric(gsub("[^0-9.]", "", timepoint_raw)),
      TRUE ~ NA_real_
    )
  )

## ---- timepoint diagnostics ----
sort(unique(samples$timepoint_raw))

samples %>%
  filter(is.na(timepoint_hr)) %>%
  count(timepoint_raw, sort = TRUE)

## ----------------
## 4. Define condition and control
## ----------------
sort(unique(samples$delivery_type))

samples <- samples %>%
  mutate(
    delivery_type = trimws(delivery_type),
    condition = case_when(
      tolower(delivery_type) == "baseline" ~ "control",
      TRUE ~ "delivery"
    ),
    control_type = if_else(condition == "control", "naive", NA_character_)
  )

# confirm condition assignment
table(samples$condition, useNA = "ifany")
table(samples$delivery_type, samples$condition, useNA = "ifany")

## ----------------
## 5. Final structural checks
## ----------------
# baseline controls
samples %>%
  filter(condition == "control") %>%
  count(tissue, timepoint_hr)

# delivery samples
samples %>%
  filter(condition == "delivery") %>%
  count(tissue, delivery_type, timepoint_hr)

#save
write_tsv(
  samples,
  file.path(
    "C:/Users/john/Desktop/IMRS_Project/00_metadata",
    paste0(dataset_id, "_samples.tsv")
  )
)
