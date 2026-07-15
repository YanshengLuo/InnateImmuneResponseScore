#!/usr/bin/env Rscript
# ============================================================
# GSE139529 — Fill missing group in design.tsv by:
#   GSM (from GEO website) -> group_label
#   GSM <-> SRR mapping via SraRunTable.csv
# Then update design.tsv "group" (or create it) using SRR (sample_id).
#
# INPUTS (edit paths):
#   - SraRunTable.csv
#   - GSE139529_design.tsv
#
# OUTPUT:
#   - GSE139529_design.filled_group.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# -------------------------
# USER SETTINGS (Windows)
# -------------------------
base_dir <- "D:/IMRS_Project/00_metadata/validations/GSE139529"
run_table_csv <- file.path(base_dir, "SraRunTable.csv")
design_tsv    <- file.path(base_dir, "GSE139529_design.tsv")
out_tsv       <- file.path(base_dir, "GSE139529_design.filled_group.tsv")

# -------------------------
# GSM -> group mapping (from GEO website)
# -------------------------
gsm_map <- tribble(
  ~GSM,        ~group_raw,
  "GSM4142998","Ad_fHbp",
  "GSM4142999","Ad_fHbp",
  "GSM4143000","Ad_fHbp",
  "GSM4143001","Ad_fHbp",
  "GSM4143002","Ad_fHbp",
  "GSM4143003","Ad_fHbp",
  "GSM4143004","MVA_fHbp",
  "GSM4143005","MVA_fHbp",
  "GSM4143006","MVA_fHbp",
  "GSM4143007","MVA_fHbp",
  "GSM4143008","MVA_fHbp",
  "GSM4143009","MVA_fHbp",
  "GSM4143010","Bexsero",
  "GSM4143011","Bexsero",
  "GSM4143012","Bexsero",
  "GSM4143013","Bexsero",
  "GSM4143014","Bexsero",
  "GSM4143015","Bexsero",
  "GSM4143016","Bex.Bex",
  "GSM4143017","Bex.Bex",
  "GSM4143018","Bex.Bex",
  "GSM4143019","Bex.Bex",
  "GSM4143020","Bex.Bex",
  "GSM4143021","Bex.Bex",
  "GSM4143022","MVA.Ad",
  "GSM4143023","MVA.Ad",
  "GSM4143024","MVA.Ad",
  "GSM4143025","MVA.Ad",
  "GSM4143026","MVA.Ad",
  "GSM4143027","MVA.Ad",
  "GSM4143028","None",
  "GSM4143029","None",
  "GSM4143030","None",
  "GSM4143031","None",
  "GSM4143032","None",
  "GSM4143033","None",
  "GSM4143034","Ad.MVA",
  "GSM4143035","Ad.MVA",
  "GSM4143036","Ad.MVA",
  "GSM4143037","Ad.MVA",
  "GSM4143038","Ad.MVA",
  "GSM4143039","Ad.MVA",
  "GSM4143040","Ad_empty",
  "GSM4143041","Ad_empty",
  "GSM4143042","Ad_empty",
  "GSM4143043","Ad_empty",
  "GSM4143044","Ad_empty",
  "GSM4143045","Ad_empty",
  "GSM4143046","Trumenba",
  "GSM4143047","Trumenba",
  "GSM4143048","Trumenba",
  "GSM4143049","Trumenba",
  "GSM4143050","Trumenba",
  "GSM4143051","Trumenba",
  "GSM4143052","Tru.Tru",
  "GSM4143053","Tru.Tru",
  "GSM4143054","Tru.Tru",
  "GSM4143055","Tru.Tru",
  "GSM4143056","Tru.Tru",
  "GSM4143057","Ad_preF",
  "GSM4143058","Ad_preF",
  "GSM4143059","Ad_preF",
  "GSM4143060","Ad_preF",
  "GSM4143061","Ad_preF",
  "GSM4143062","Ad_preF"
)

sanitize_group <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("[/]", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\.]+", "_") %>%
    str_replace_all("\\.+", ".") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

infer_condition_simple <- function(group_raw) {
  # Rule: "None" is CONTROL, everything else DELIVERY
  ifelse(str_to_lower(group_raw) == "none", "CONTROL", "DELIVERY")
}

# -------------------------
# Load files
# -------------------------
if (!file.exists(run_table_csv)) stop("Missing: ", run_table_csv)
if (!file.exists(design_tsv))    stop("Missing: ", design_tsv)

rt <- read_csv(run_table_csv, show_col_types = FALSE)
des <- read_tsv(design_tsv, show_col_types = FALSE)

# -------------------------
# Identify columns in SraRunTable
# -------------------------
# SRR column
run_col_candidates <- c("Run", "run", "SRR", "srr", "srr_id")
run_col <- run_col_candidates[run_col_candidates %in% colnames(rt)][1]
if (is.na(run_col)) {
  stop("Could not find SRR/Run column in SraRunTable.csv. Columns are: ",
       paste(colnames(rt), collapse = ", "))
}

# GSM column (GEO accession)
gsm_col_candidates <- c("GEO_Accession", "geo_accession", "GSM", "gsm", "GEO", "Sample_Name", "sample_name")
gsm_col <- gsm_col_candidates[gsm_col_candidates %in% colnames(rt)][1]
if (is.na(gsm_col)) {
  # Try fuzzy find: anything containing "GEO" or "GSM"
  idx <- which(str_detect(colnames(rt), regex("geo|gsm", ignore_case = TRUE)))
  if (length(idx) == 1) {
    gsm_col <- colnames(rt)[idx]
  } else {
    stop("Could not confidently find GSM/GEO column in SraRunTable.csv. Columns are: ",
         paste(colnames(rt), collapse = ", "),
         "\nTip: Open SraRunTable.csv and tell me the exact column name containing GSM IDs.")
  }
}

# -------------------------
# Build SRR -> GSM mapping
# -------------------------
srr2gsm <- rt %>%
  transmute(
    sample_id = as.character(.data[[run_col]]),
    GSM       = as.character(.data[[gsm_col]])
  ) %>%
  filter(!is.na(sample_id), sample_id != "", !is.na(GSM), GSM != "") %>%
  distinct()

# -------------------------
# Join GSM -> group, then update design
# -------------------------
mapped <- srr2gsm %>%
  left_join(gsm_map, by = "GSM") %>%
  mutate(
    group = sanitize_group(group_raw),
    condition_simple = infer_condition_simple(group_raw)
  )

# sanity: any SRR with missing group?
missing_group <- mapped %>% filter(is.na(group_raw) | group_raw == "")
if (nrow(missing_group) > 0) {
  stop("Some SRRs could not be mapped to a group via GSM. Example rows:\n",
       paste(utils::capture.output(print(head(missing_group, 10))), collapse = "\n"))
}

# ensure design has sample_id
if (!("sample_id" %in% colnames(des))) {
  stop("design.tsv must contain a 'sample_id' column (SRR IDs). Found columns: ",
       paste(colnames(des), collapse = ", "))
}

# merge into design
des_out <- des %>%
  left_join(mapped %>% select(sample_id, group, condition_simple), by = "sample_id")

# if design already has 'group'/'condition_simple', replace them
if ("group.x" %in% colnames(des_out) && "group.y" %in% colnames(des_out)) {
  des_out <- des_out %>%
    mutate(group = ifelse(!is.na(group.y), group.y, group.x)) %>%
    select(-group.x, -group.y)
} else if ("group" %in% colnames(des_out)) {
  # keep as is (already named group)
} else {
  # join created group column directly
}

if ("condition_simple.x" %in% colnames(des_out) && "condition_simple.y" %in% colnames(des_out)) {
  des_out <- des_out %>%
    mutate(condition_simple = ifelse(!is.na(condition_simple.y), condition_simple.y, condition_simple.x)) %>%
    select(-condition_simple.x, -condition_simple.y)
}

# hard check: all sample_ids in design mapped?
unmapped_design <- des_out %>% filter(is.na(group) | group == "")
if (nrow(unmapped_design) > 0) {
  stop("Some sample_id values in design.tsv did not map to GSM/group. Example:\n",
       paste(utils::capture.output(print(head(unmapped_design, 10))), collapse = "\n"))
}

# -------------------------
# Write
# -------------------------
write_tsv(des_out, out_tsv)
message("Wrote: ", out_tsv)

# Optional: print counts
message("\nGroup counts:")
print(des_out %>% count(group, sort = TRUE), n = 50)

message("\nCondition_simple counts:")
print(des_out %>% count(condition_simple, sort = TRUE))