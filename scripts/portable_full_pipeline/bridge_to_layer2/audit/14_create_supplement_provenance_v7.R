suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

root <- normalizePath(Sys.getenv("IMRS_PORTED_OUTPUT_ROOT",
                                 Sys.getenv("IMRS_PROJECT_ROOT", ".")),
                      winslash = "/", mustWork = FALSE)
verified_dir <- normalizePath(Sys.getenv("IMRS_BRIDGE_VERIFIED_METADATA_ROOT",
                                         file.path(root, "00_metadata", "verified_metadata")),
                              winslash = "/", mustWork = FALSE)
results_dir <- normalizePath(Sys.getenv("IMRS_LAYER2_AUDIT_RESULTS_DIR",
                                        file.path(root, "audit", "results")),
                             winslash = "/", mustWork = FALSE)
report_dir <- normalizePath(Sys.getenv("IMRS_LAYER2_AUDIT_REPORT_DIR",
                                       file.path(root, "audit", "report")),
                            winslash = "/", mustWork = FALSE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  role_table = file.path(results_dir, "manuscript_dataset_role_table.tsv"),
  weak_context = file.path(results_dir, "weak_dataset_paper_context_audit.tsv"),
  weak_interp = file.path(results_dir, "manuscript_weak_dataset_interpretation_table.tsv"),
  dataset_audit = file.path(results_dir, "dataset_audit_table.tsv"),
  step09_eval = file.path(root, "05_score", "transfer", "eval", "step09_split_eval.tsv"),
  step09_summary = file.path(root, "05_score", "transfer", "eval", "step09_split_summary.tsv"),
  step09_samples = file.path(root, "05_score", "transfer", "eval", "step09_split_sample_level.tsv")
)

out_paths <- list(
  folder_inventory = file.path(results_dir, "all_verified_metadata_folder_inventory.tsv"),
  split_inventory = file.path(results_dir, "all_split_design_inventory.tsv"),
  tsv = file.path(results_dir, "supplement_dataset_split_provenance_v7.tsv"),
  csv = file.path(results_dir, "supplement_dataset_split_provenance_v7.csv"),
  xlsx = file.path(results_dir, "supplement_dataset_split_provenance_v7.xlsx"),
  rmd = file.path(report_dir, "Supplement_Dataset_Split_Provenance_v7.Rmd"),
  pdf = file.path(report_dir, "Supplement_Dataset_Split_Provenance_v7.pdf")
)

expected_metadata_folders <- c(
  "GSE39129", "GSE119119", "GSE139529", "GSE166655", "GSE167521",
  "GSE178313", "GSE190850_HUMAN", "GSE237068", "GSE262515",
  "GSE264344", "GSE279372", "GSE279743", "GSE279744", "GSE314070"
)

read_tsv_safe <- function(path) {
  if (!file.exists(path)) {
    warning("Missing input: ", path)
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "N/A"))
}

read_any_table <- function(path) {
  if (!file.exists(path)) return(tibble())
  if (str_detect(tolower(path), "\\.csv$")) {
    suppressWarnings(read_csv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "N/A")))
  } else {
    suppressWarnings(read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "N/A")))
  }
}

as_chr <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "N/A")] <- NA_character_
  x
}

coalesce_chr <- function(...) {
  vals <- lapply(list(...), as_chr)
  if (length(vals) == 0) return(character())
  out <- vals[[1]]
  if (length(vals) > 1) {
    for (i in seq(2, length(vals))) out <- dplyr::coalesce(out, vals[[i]])
  }
  as_chr(out)
}

collapse_unique <- function(x, max_items = Inf) {
  x <- unique(as_chr(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  if (is.finite(max_items) && length(x) > max_items) {
    return(paste0(paste(head(x, max_items), collapse = "; "), "; ..."))
  }
  paste(x, collapse = "; ")
}

pick_col <- function(df, candidates) {
  if (length(names(df)) == 0) return(NA_character_)
  lower_names <- tolower(names(df))
  for (candidate in candidates) {
    hit <- which(lower_names == tolower(candidate))
    if (length(hit) > 0) return(names(df)[hit[[1]]])
  }
  for (candidate in candidates) {
    hit <- which(str_detect(lower_names, fixed(tolower(candidate))))
    if (length(hit) > 0) return(names(df)[hit[[1]]])
  }
  NA_character_
}

ensure_columns <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df
}

truthy <- function(x) {
  x <- tolower(as.character(x))
  x %in% c("true", "t", "1", "yes", "y")
}

clean_label <- function(x) {
  out <- as_chr(x)
  out <- str_replace_all(out, "_", " ")
  str_squish(out)
}

short_text <- function(x, n = 80) {
  x <- as_chr(x)
  str_trunc(x, n)
}

base_gse <- function(x) str_extract(as_chr(x), "GSE[0-9]+")

pub_key_from_dataset <- function(dataset_id, gse_id = NA_character_) {
  out <- coalesce_chr(base_gse(dataset_id), base_gse(gse_id))
  out
}

parse_token <- function(x, key) {
  x <- tools::file_path_sans_ext(basename(as_chr(x)))
  pattern <- switch(
    key,
    T = "__T=(.*?)__H=",
    H = "__H=(.*?)__B=",
    B = "__B=(.*?)__G=",
    G = "__G=(.*?)__VS=",
    VS = "__VS=(.*)$"
  )
  as_chr(str_match(x, pattern)[, 2])
}

parse_time <- function(split_id, split_file = NA_character_) {
  h <- parse_token(split_id, "H")
  out <- suppressWarnings(as.numeric(h))
  idx <- is.na(out)
  if (any(idx)) {
    src <- paste(as_chr(split_id), basename(as_chr(split_file)), sep = " ")
    m <- str_match(src, regex("(?:__H=|\\bH=|time[_-]?)([0-9.]+)|([0-9.]+)h", ignore_case = TRUE))
    parsed <- suppressWarnings(as.numeric(coalesce(m[, 2], m[, 3])))
    out[idx] <- parsed[idx]
  }
  out
}

parse_tissue <- function(split_id, split_file = NA_character_) {
  out <- parse_token(split_id, "T")
  idx <- is.na(out)
  if (any(idx)) {
    src <- paste(as_chr(split_id), basename(as_chr(split_file)), sep = " ")
    parsed <- str_match(src, regex("(?:__T=|tissue[_-])([^_]+)", ignore_case = TRUE))[, 2]
    out[idx] <- parsed[idx]
  }
  clean_label(out)
}

role_labels <- c(
  anchor = "Locked anchor",
  anchor_discovery = "Locked anchor",
  strict_platform_anchor = "Locked anchor",
  strict_anchor = "Locked anchor",
  additional_acute_discovery_support = "Locked anchor",
  discovery_support = "Locked anchor",
  external_acute = "Primary acute validation",
  primary_acute_validation = "Primary acute validation",
  external_extended = "Extended validation",
  extended_validation = "Extended validation",
  calibration = "Secondary support",
  secondary_support = "Secondary support",
  secondary_support_not_primary = "Secondary support",
  excluded_or_unclear = "Excluded/unclear"
)

manual_paper_context_mapping <- tribble(
  ~dataset_id_or_gse_id, ~dataset_folder_alias, ~PMID_manual, ~DOI_manual, ~publication_year_manual, ~paper_title_manual, ~organism_manual, ~source_database_manual, ~accession_manual, ~delivery_platform_manual, ~paper_treatment_context, ~expected_imrs_behavior, ~paper_result_alignment_rule, ~manuscript_context_note, ~limitations_note,
  "GSE39129", NA, "22871668", "10.1038/mt.2012.150", 2012L, "A TLR and Non-TLR Mediated Innate Response to Lentiviruses Restricts Hepatocyte Entry and Can be Ameliorated by Pharmacological Blockade", "Mus musculus", "GEO", "GSE39129", "lentiviral vector", "Mouse liver RNA-seq 4 h after systemic lentiviral vector administration versus saline/control. Paper reports liver innate and interferon response after LV administration.", "strong_expected_acute_innate_signal", "Strong positive IMRS is expected because this is an acute 4 h liver response to lentiviral vector delivery.", "Locked anchor / acute viral-vector innate response. Good anchor evidence.", "Mouse liver vector response; not clinical reactogenicity.",
  "GSE119119", NA, "30209217", "10.1073/pnas.1806314115", 2018L, "TRIM21 mediates antibody inhibition of adenovirus-based gene delivery and vaccination", "Mus musculus", "GEO", "GSE119119", "adenoviral vector / Ad5 antibody perturbation", "Mouse liver RNA-seq 4 h after Ad5 delivery and Ad5 plus antibody perturbation in WT and TRIM21 KO contexts. Current splits compare virus/antibody perturbation groups against WT uninfected whole liver.", "strong_expected_acute_innate_signal", "Positive IMRS is expected because this is an acute adenoviral-vector perturbation with innate antiviral biology at 4 h. Genotype/antibody perturbation should be mentioned as context.", "Primary acute validation with perturbation context.", "Not a clean formulation-only comparison because TRIM21 genotype and antibody perturbation are involved.",
  "GSE139529", NA, "34717548", "10.1186/s12864-021-08061-8", 2021L, "Distinct patterns of whole blood transcriptional responses are induced in mice following immunisation with adenoviral and poxviral vector vaccines encoding the same antigen", "Mus musculus", "GEO", "GSE139529", "adenoviral and modified vaccinia Ankara viral vector vaccines", "Mouse whole-blood RNA-seq 24 h after immunisation with adenoviral and poxviral vector vaccines encoding the same antigen.", "strong_expected_acute_innate_signal", "Strong positive IMRS is expected for most vector vaccine groups because the paper studies acute blood transcriptional responses after viral-vector immunisation.", "Primary acute validation.", "Whole blood response includes systemic immune-cell composition and vector-specific antigen effects.",
  "GSE166655", NA, "34417184", "10.1126/sciadv.abi6896", 2021L, "Gene therapy with AR isoform 2 rescues spinal and bulbar muscular atrophy phenotype by modulating AR transcriptional activity", "Mus musculus", "GEO", "GSE166655", "AAV9 gene therapy", "AAV9-AR45 or AAV9-GFP gene therapy in SBMA mouse model; RNA-seq from tibialis anterior muscle 1008 h (6 weeks) after AAV9 injection. Current splits are AAV9-AR45 or AAV9-GFP versus untreated.", "weak_but_biologically_explainable", "Weak or near-zero acute IMRS is expected because this is 1008 h after delivery and dominated by therapeutic transgene/disease-rescue biology rather than acute innate activation.", "Extended validation only, not primary acute validation.", "Late AAV therapeutic context; not suitable as strict acute innate validation.",
  "GSE167521", NA, "34841223", "10.1016/j.isci.2021.103479", 2021L, "The mRNA-LNP platform's lipid nanoparticle component used in preclinical vaccine studies is highly inflammatory", "Mus musculus", "GEO", "GSE167521", "lipid nanoparticle", "Mouse LNP exposure dataset. The paper reports that preclinical LNPs induce strong inflammatory responses, including neutrophil infiltration, inflammatory pathways, cytokines, and chemokines.", "strong_expected_acute_innate_signal", "Positive IMRS is expected because the study directly describes inflammatory innate activation from LNP exposure.", "Locked anchor / acute LNP inflammatory response.", "Preclinical LNP inflammation model; not a clinical safety endpoint.",
  "GSE178313", NA, "35132167", "10.1038/s41565-021-01030-y", 2022L, "Species-dependent in vivo mRNA delivery and cellular responses to nanoparticles", "Homo sapiens; Mus musculus", "GEO", "GSE178313", "lipid nanoparticle", "Species-dependent LNP delivery and cellular-response study; transcriptomic hepatocyte response was measured 24 h after LNP or PBS injection in humanized/mouse liver context.", "moderate_expected_signal", "Moderate positive IMRS is biologically plausible because the paper studies 24-h LNP-induced cellular responses; interpretation remains extended/context-shifted because of the species/humanized context.", "Extended validation, not primary acute validation.", "Species-transfer context; tissue/cell composition may differ from mouse acute anchors.",
  "GSE190850", "GSE190850_HUMAN", "34990408", "10.1172/jci.insight.155655", 2022L, "Aerosol delivery, but not intramuscular injection, of adenovirus-vectored tuberculosis vaccine induces respiratory-mucosal immunity in humans", "Homo sapiens", "GEO", "GSE190850", "human adenovirus-vectored tuberculosis vaccine", "Human phase Ib AdHu5Ag85A adenovirus-vectored tuberculosis vaccine study comparing inhaled aerosol and intramuscular delivery. Samples include peripheral blood and bronchoalveolar lavage immune readouts.", "metadata_only_not_interpreted", "If not scored in current IMRS evaluation, keep as metadata-only. If scored later, interpret cautiously as human vaccine route/mucosal-immunity context, not mouse acute anchor behavior.", "Metadata-only / possible future human validation.", "Human clinical design, route-specific mucosal immunity, may not match mouse delivery-induced innate anchors.",
  "GSE237068", NA, "38101397", "10.1016/j.cmet.2023.11.010", 2024L, "ChREBP is activated by reductive stress and mediates GCKR-associated metabolic traits", "Mus musculus", "GEO", "GSE237068", "adenoviral vector (EcSTH or GFP)", "Metabolic/reductive-stress study using retro-orbital EcSTH-expressing or GFP-control adenovirus; liver RNA-seq was collected 4 days (96 h) after administration.", "excluded_or_not_delivery_relevant", "Keep as metadata-only/excluded because the study is a metabolic EcSTH perturbation rather than a clean delivery-induced innate activation experiment.", "Metadata-audited but excluded/contextual.", "Verified liver adenoviral-vector metadata; no current IMRS score, so no delivery-response claim is made.",
  "GSE262515", NA, "38902241", "10.1038/s41392-024-01871-8", 2024L, "Optimized RNA interference therapeutics combined with interleukin-2 mRNA for treating hepatitis B virus infection", "Homo sapiens; Mus musculus", "GEO", "GSE262515", "optimized therapeutic lipid nanoparticle / tLNP RNAi", "Therapeutic tLNP/RNAi context with a 16-h HepG2 cell-line arm and a mouse tissue/liver arm recorded as 72 h by GEO and supplementary RNA-seq Methods. The main-paper Figure 4e caption instead says 48 h for the mouse arm, so that arm remains explicitly source-conflicted.", "weak_but_biologically_explainable", "Weak or variable IMRS is explainable by therapeutic, cargo/formulation-specific, and often baseline-vs-treatment contrasts rather than a clean shared-vehicle acute delivery contrast; the cell-line arm is not a late timepoint.", "Secondary support / non-primary evidence.", "Do not claim direct shared-vehicle comparison unless the split is explicitly tLNP/siHBV versus tLNP/siNC; retain the mouse-arm 72-h/48-h source conflict disclosure.",
  "GSE264344", NA, "40162786", "10.1128/jvi.00247-25", 2025L, "Early spatiotemporal evolution of the immune response elicited by adenovirus serotype 26 vector vaccination in mice", "Mus musculus", "GEO", "GSE264344", "adenoviral vector", "Mouse Ad26/Ad5 vector vaccination time course across blood, muscle, and draining lymph node from 1 h to 72 h. The paper reports that innate immune response begins by 1 h, evolves across tissues within 24 h, and wanes by 72 h.", "strong_expected_acute_innate_signal_for_1_to_24h; weak_but_biologically_explainable_for_72h", "Strong positive IMRS at 1-24 h matches acute innate tissue response. Weak 72 h IMRS matches waning tissue-time kinetics, not model failure.", "1-24 h rows = locked anchor; 72 h rows = extended validation/kinetic context.", "Time-course rows must not be pooled as one biological role.",
  "GSE279372", NA, "39494905", "10.1128/jvi.01699-24", 2024L, "Lipid nanoparticles as adjuvant of norovirus VLP vaccine augment cellular and humoral immune responses in a TLR9- and type I IFN-dependent pathway", "Mus musculus", "GEO", "GSE279372", "LNP adjuvant / norovirus VLP vaccine", "Mouse draining lymph node RNA-seq 24 h after NoV VLP vaccine and LNP/adjuvant treatment. Paper reports LNP adjuvant activity and type I IFN/TLR9-dependent mechanism.", "strong_expected_acute_innate_signal", "Strong positive IMRS is biologically plausible because the LNP adjuvant context is associated with innate/type I IFN-dependent vaccine response.", "Locked anchor evidence in the five-dataset production anchor set.", "VLP antigen and adjuvant biology are mixed.",
  "GSE279743", NA, "39690742", "10.1016/j.ymthe.2024.12.032", 2025L, "Low-inflammatory lipid nanoparticle-based mRNA vaccine elicits protective immunity against H5N1 influenza virus with reduced adverse reactions", "Mus musculus", "GEO", "GSE279743", "mRNA lipid nanoparticle vaccine", "Whole lymph node RNA-seq after PBS, LNPssPalmO, and/or LNPSM-102 treatment in a low-inflammatory H5N1 mRNA-LNP vaccine study.", "moderate_expected_signal", "Positive IMRS is expected for acute LNP vaccine contrasts, but lower signals for low-inflammatory formulations are biologically plausible because the paper explicitly focuses on reduced inflammatory/adverse-reaction LNP design.", "Primary acute validation.", "Formulation intentionally reduces inflammation; contrast-specific interpretation matters.",
  "GSE279744", NA, "39690742", "10.1016/j.ymthe.2024.12.032", 2025L, "Low-inflammatory lipid nanoparticle-based mRNA vaccine elicits protective immunity against H5N1 influenza virus with reduced adverse reactions", "Mus musculus", "GEO", "GSE279744", "mRNA lipid nanoparticle vaccine", "Lymph node dendritic-cell RNA-seq from the same H5N1 low-inflammatory LNP vaccine study.", "strong_expected_acute_innate_signal", "Strong IMRS in sorted lymph-node dendritic-cell context is plausible because dendritic cells are innate immune responders, and these acute rows are included in the production model.", "Locked anchor evidence in the five-dataset production anchor set.", "Cell-sorted immune population from same study; interpret separately from primary validation.",
  "GSE314070", NA, "41671339", "10.1126/scitranslmed.adw6105", 2026L, "Parenteral vaccination with an adjuvanted mRNA vaccine induces protective mucosal immunity against rotavirus in neonatal mice", "Mus musculus", "GEO", "GSE314070", "adjuvanted mRNA lipid nanoparticle vaccine / Am80-LNP", "Small intestine RNA-seq at 336 h/day 14 after booster immunization with Am80-LNP rotavirus mRNA vaccine. The platform is designed to induce gut mucosal immunity and mitigate injection-site inflammation.", "weak_but_biologically_explainable", "Weak or variable acute IMRS is expected because this is late small-intestine mucosal/adaptive biology after booster, not acute injection-site innate activation.", "Extended validation or supplement.", "Late distal-tissue adaptive/mucosal context; not primary acute validation."
) %>%
  mutate(pub_key = dataset_id_or_gse_id)

role <- read_tsv_safe(paths$role_table)
weak_context <- read_tsv_safe(paths$weak_context)
weak_interp <- read_tsv_safe(paths$weak_interp)
dataset_audit <- read_tsv_safe(paths$dataset_audit)
eval <- read_tsv_safe(paths$step09_eval)
sample_level <- read_tsv_safe(paths$step09_samples)

if (nrow(role) == 0) stop("manuscript_dataset_role_table.tsv is required.")

metadata_files <- list.files(
  verified_dir,
  pattern = "\\.(tsv|csv)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
) %>%
  tibble(path = .) %>%
  mutate(
    path_norm = str_replace_all(path, "\\\\", "/"),
    file_name = basename(path_norm),
    lower_path = tolower(path_norm),
    lower_file = tolower(file_name),
    under_split_dir = str_detect(lower_path, "/split(?:t)?ed/"),
    under_scoring_dir = str_detect(lower_path, "/scoring/"),
    is_index_or_map = str_detect(lower_file, "contrast_index|group_map|control_report|control_audit"),
    is_split_file = under_split_dir & !is_index_or_map & str_detect(file_name, "__G=.*__VS=.*\\.(tsv|csv)$"),
    is_design_file = !is_split_file & str_detect(file_name, "_design\\.(tsv|csv)$"),
    dataset_key_from_file = str_remove(file_name, "_design\\.(tsv|csv)$"),
    dataset_key_base = base_gse(dataset_key_from_file)
  )

design_files <- metadata_files %>%
  filter(is_design_file) %>%
  mutate(
    design_key = str_remove(file_name, "_design\\.(tsv|csv)$"),
    priority = case_when(
      under_scoring_dir ~ 1L,
      TRUE ~ 2L
    )
  )

find_design_file <- function(dataset_id, gse_id = NA_character_) {
  keys <- unique(na.omit(c(
    as_chr(dataset_id),
    as_chr(gse_id),
    base_gse(dataset_id),
    base_gse(gse_id)
  )))
  candidates <- design_files %>%
    filter(design_key %in% keys | dataset_key_base %in% keys) %>%
    arrange(priority, nchar(path_norm))
  if (nrow(candidates) == 0) return(NA_character_)
  candidates$path_norm[[1]]
}

split_file_index <- metadata_files %>%
  filter(is_split_file) %>%
  transmute(
    source_split_file = path_norm,
    split_id = tools::file_path_sans_ext(file_name),
    split_base = split_id,
    folder = basename(dirname(path_norm)),
    gse_id = coalesce_chr(base_gse(folder), base_gse(split_id)),
    dataset_id_guess = case_when(
      str_detect(folder, "GSE262515_cell_line") ~ "GSE262515_cell_line",
      str_detect(folder, "GSE262515_tissue") ~ "GSE262515_tissue",
      str_detect(folder, "GSE190850_HUMAN") ~ "GSE190850_HUMAN",
      TRUE ~ coalesce_chr(str_remove(folder, "_design$"), gse_id)
    ),
    tissue_from_filename = parse_tissue(split_id, source_split_file),
    time_h_from_filename = parse_time(split_id, source_split_file),
    treatment_token = parse_token(split_id, "G"),
    control_token = parse_token(split_id, "VS")
  )

summarise_split_file <- function(split_path, split_id) {
  if (is.na(split_path) || !file.exists(split_path)) {
    return(tibble(
      source_split_file = as_chr(split_path),
      sample_ids_delivery = NA_character_,
      sample_ids_control = NA_character_,
      n_delivery_meta = NA_integer_,
      n_control_meta = NA_integer_,
      condition_simple_meta = NA_character_,
      treatment_group_raw_meta = NA_character_,
      control_group_raw_meta = NA_character_,
      tissue_meta = parse_tissue(split_id, split_path),
      time_h_meta = parse_time(split_id, split_path),
      timepoint_original = NA_character_,
      route_meta = NA_character_,
      dose_meta = NA_character_,
      platform_meta = NA_character_,
      vector_meta = NA_character_,
      formulation_meta = NA_character_,
      cargo_meta = NA_character_,
      transgene_meta = NA_character_,
      organism_meta = NA_character_,
      cell_line_meta = NA_character_
    ))
  }

  df <- read_any_table(split_path)
  if (nrow(df) == 0) {
    return(tibble(
      source_split_file = str_replace_all(split_path, "\\\\", "/"),
      sample_ids_delivery = NA_character_,
      sample_ids_control = NA_character_,
      n_delivery_meta = 0L,
      n_control_meta = 0L,
      condition_simple_meta = NA_character_,
      treatment_group_raw_meta = NA_character_,
      control_group_raw_meta = NA_character_,
      tissue_meta = parse_tissue(split_id, split_path),
      time_h_meta = parse_time(split_id, split_path),
      timepoint_original = NA_character_,
      route_meta = NA_character_,
      dose_meta = NA_character_,
      platform_meta = NA_character_,
      vector_meta = NA_character_,
      formulation_meta = NA_character_,
      cargo_meta = NA_character_,
      transgene_meta = NA_character_,
      organism_meta = NA_character_,
      cell_line_meta = NA_character_
    ))
  }

  sample_col <- pick_col(df, c("sample_id", "sample", "run", "geo_accession"))
  condition_col <- pick_col(df, c("condition_simple", "condition", "class"))
  is_control_col <- pick_col(df, c("is_control"))
  group_raw_col <- pick_col(df, c("group_raw", "group", "treatment", "condition_name"))
  treatment_col <- pick_col(df, c("treatment", "treatment_group"))
  control_col <- pick_col(df, c("control", "control_group", "control_label"))
  tissue_col <- pick_col(df, c("tissue_raw", "tissue", "organ"))
  cell_line_col <- pick_col(df, c("cell_line", "cell_type", "celltype", "cell"))
  organism_col <- pick_col(df, c("organism", "species"))
  time_col <- pick_col(df, c("time_h", "timepoint_hr", "timepoint", "time_raw", "time"))
  route_col <- pick_col(df, c("route"))
  dose_col <- pick_col(df, c("dose"))
  platform_col <- pick_col(df, c("platform", "delivery_platform"))
  vector_col <- pick_col(df, c("vector"))
  formulation_col <- pick_col(df, c("formulation"))
  cargo_col <- pick_col(df, c("cargo"))
  transgene_col <- pick_col(df, c("transgene"))

  condition <- if (!is.na(condition_col)) as_chr(df[[condition_col]]) else rep(NA_character_, nrow(df))
  is_control <- if (!is.na(is_control_col)) truthy(df[[is_control_col]]) else rep(NA, nrow(df))
  group_raw <- if (!is.na(group_raw_col)) as_chr(df[[group_raw_col]]) else rep(NA_character_, nrow(df))
  control_from_split <- parse_token(split_id, "VS")
  group_from_split <- parse_token(split_id, "G")

  control_rows <- ifelse(!is.na(condition), toupper(condition) == "CONTROL", NA)
  control_rows <- ifelse(is.na(control_rows) & !is.na(is_control), is_control, control_rows)
  control_rows <- ifelse(
    is.na(control_rows) & !is.na(group_raw) & !is.na(control_from_split),
    group_raw == control_from_split,
    control_rows
  )
  control_rows[is.na(control_rows)] <- FALSE

  delivery_rows <- ifelse(!is.na(condition), toupper(condition) == "DELIVERY", NA)
  delivery_rows <- ifelse(is.na(delivery_rows) & !is.na(is_control), !is_control, delivery_rows)
  delivery_rows <- ifelse(
    is.na(delivery_rows) & !is.na(group_raw) & !is.na(group_from_split),
    group_raw == group_from_split,
    delivery_rows
  )
  delivery_rows[is.na(delivery_rows)] <- FALSE

  sample_ids <- if (!is.na(sample_col)) as_chr(df[[sample_col]]) else rep(NA_character_, nrow(df))
  time_vals <- if (!is.na(time_col)) as_chr(df[[time_col]]) else rep(NA_character_, nrow(df))
  active_rows <- delivery_rows | control_rows
  if (!any(active_rows)) active_rows <- rep(TRUE, nrow(df))

  tibble(
    source_split_file = str_replace_all(split_path, "\\\\", "/"),
    sample_ids_delivery = collapse_unique(sample_ids[delivery_rows]),
    sample_ids_control = collapse_unique(sample_ids[control_rows]),
    n_delivery_meta = sum(delivery_rows, na.rm = TRUE),
    n_control_meta = sum(control_rows, na.rm = TRUE),
    condition_simple_meta = collapse_unique(condition),
    treatment_group_raw_meta = collapse_unique(coalesce_chr(
      if (!is.na(treatment_col)) df[[treatment_col]][delivery_rows] else NA_character_,
      group_raw[delivery_rows],
      group_from_split
    )),
    control_group_raw_meta = collapse_unique(coalesce_chr(
      if (!is.na(control_col)) df[[control_col]][control_rows] else NA_character_,
      group_raw[control_rows],
      control_from_split
    )),
    tissue_meta = collapse_unique(if (!is.na(tissue_col)) df[[tissue_col]][active_rows] else parse_tissue(split_id, split_path)),
    time_h_meta = suppressWarnings(as.numeric(collapse_unique(if (!is.na(time_col)) time_vals[active_rows] else parse_time(split_id, split_path), max_items = 1))),
    timepoint_original = collapse_unique(time_vals[active_rows]),
    route_meta = collapse_unique(if (!is.na(route_col)) df[[route_col]][active_rows] else NA_character_),
    dose_meta = collapse_unique(if (!is.na(dose_col)) df[[dose_col]][active_rows] else NA_character_),
    platform_meta = collapse_unique(if (!is.na(platform_col)) df[[platform_col]][active_rows] else NA_character_),
    vector_meta = collapse_unique(if (!is.na(vector_col)) df[[vector_col]][active_rows] else NA_character_),
    formulation_meta = collapse_unique(if (!is.na(formulation_col)) df[[formulation_col]][active_rows] else NA_character_),
    cargo_meta = collapse_unique(if (!is.na(cargo_col)) df[[cargo_col]][active_rows] else NA_character_),
    transgene_meta = collapse_unique(if (!is.na(transgene_col)) df[[transgene_col]][active_rows] else NA_character_),
    organism_meta = collapse_unique(if (!is.na(organism_col)) df[[organism_col]][active_rows] else NA_character_),
    cell_line_meta = collapse_unique(if (!is.na(cell_line_col)) df[[cell_line_col]][active_rows] else NA_character_)
  )
}

split_inventory <- bind_rows(lapply(seq_len(nrow(split_file_index)), function(i) {
  sm <- summarise_split_file(split_file_index$source_split_file[[i]], split_file_index$split_id[[i]])
  tibble(
    source_split_file = split_file_index$source_split_file[[i]],
    split_id = split_file_index$split_id[[i]],
    dataset_id_guess = split_file_index$dataset_id_guess[[i]],
    gse_id = split_file_index$gse_id[[i]],
    tissue_from_filename = split_file_index$tissue_from_filename[[i]],
    time_h_from_filename = split_file_index$time_h_from_filename[[i]],
    treatment_token = split_file_index$treatment_token[[i]],
    control_token = split_file_index$control_token[[i]],
    n_rows = nrow(read_any_table(split_file_index$source_split_file[[i]])),
    n_delivery = sm$n_delivery_meta[[1]],
    n_control = sm$n_control_meta[[1]],
    treatment_group_raw = sm$treatment_group_raw_meta[[1]],
    control_group_raw = sm$control_group_raw_meta[[1]],
    tissue_recovered = sm$tissue_meta[[1]],
    time_h_recovered = sm$time_h_meta[[1]],
    is_scored_in_current_eval = split_file_index$split_id[[i]] %in% role$split_id
  )
})) %>%
  arrange(dataset_id_guess, time_h_recovered, tissue_recovered, split_id)

write_tsv(split_inventory, out_paths$split_inventory, na = "NA")

role_split_paths <- role %>%
  transmute(
    split_id,
    source_split_file_role = str_replace_all(split_path, "\\\\", "/")
  ) %>%
  left_join(split_file_index %>% select(split_id, source_split_file), by = "split_id") %>%
  mutate(source_split_file_final = coalesce_chr(source_split_file_role, source_split_file))

split_meta <- bind_rows(lapply(seq_len(nrow(role_split_paths)), function(i) {
  summarise_split_file(role_split_paths$source_split_file_final[[i]], role_split_paths$split_id[[i]]) %>%
    mutate(split_id = role_split_paths$split_id[[i]], .before = 1)
}))

eval_keep <- eval %>%
  select(any_of(c(
    "split_id", "contrast_label", "control_label", "n_split_samples", "n_controls",
    "n_delivery", "delta_mean_imrs_z", "auc_imrs_z", "pass", "fail_reason",
    "delivery_sample_ids", "control_sample_ids"
  ))) %>%
  distinct(split_id, .keep_all = TRUE)

weak_keep <- weak_context %>%
  select(any_of(c(
    "dataset_id", "gse_id", "split_id", "PMID", "paper_title", "DOI",
    "publication_year", "organism", "tissue", "time_h", "treatment_group",
    "control_group", "delivery_platform", "cargo_or_transgene",
    "original_IMRS_delta", "original_AUC_secondary", "explanation_category",
    "manuscript_ready_interpretation", "reviewer_risk_level", "reviewer_risk", "notes"
  ))) %>%
  ensure_columns(c(
    "dataset_id", "gse_id", "split_id", "PMID", "paper_title", "DOI",
    "publication_year", "organism", "tissue", "time_h", "treatment_group",
    "control_group", "delivery_platform", "cargo_or_transgene",
    "original_IMRS_delta", "original_AUC_secondary", "explanation_category",
    "manuscript_ready_interpretation", "reviewer_risk_level", "reviewer_risk", "notes"
  )) %>%
  distinct(split_id, .keep_all = TRUE)

weak_interp_keep <- weak_interp %>%
  select(any_of(c(
    "split_id", "manuscript_ready_interpretation", "reviewer_risk_level",
    "reviewer_risk", "notes", "explanation_category"
  ))) %>%
  ensure_columns(c(
    "split_id", "manuscript_ready_interpretation", "reviewer_risk_level",
    "reviewer_risk", "notes", "explanation_category"
  )) %>%
  distinct(split_id, .keep_all = TRUE)

audit_keep <- dataset_audit %>%
  select(any_of(c(
    "split_id", "organism", "delivery_platform_type", "control_definition",
    "acute_late_group", "inclusion_rationale", "limitation"
  ))) %>%
  ensure_columns(c(
    "split_id", "organism", "delivery_platform_type", "control_definition",
    "acute_late_group", "inclusion_rationale", "limitation"
  )) %>%
  distinct(split_id, .keep_all = TRUE)

classify_interpretation <- function(dataset_id, manuscript_group, delta, time_h, expected) {
  dataset_id <- as_chr(dataset_id)
  manuscript_group <- as_chr(manuscript_group)
  expected <- as_chr(expected)
  case_when(
    is.na(delta) & expected == "excluded_or_not_delivery_relevant" ~ "excluded_or_not_delivery_relevant",
    is.na(delta) ~ "metadata_only_not_interpreted",
    str_detect(dataset_id, "GSE166655|GSE262515|GSE314070") ~ "weak_but_biologically_explainable",
    str_detect(dataset_id, "GSE264344") & !is.na(time_h) & time_h == 72 ~ "weak_but_biologically_explainable",
    delta <= 0 ~ "negative_or_low_requires_review",
    manuscript_group == "Primary acute validation" ~ "strong_expected_acute_innate_signal",
    str_detect(dataset_id, "GSE264344") & !is.na(time_h) & time_h <= 24 & delta > 0 ~ "strong_expected_acute_innate_signal",
    expected == "moderate_expected_signal" ~ "moderate_expected_signal",
    expected == "excluded_or_not_delivery_relevant" ~ "excluded_or_not_delivery_relevant",
    delta > 0 & delta < 1 ~ "moderate_expected_signal",
    TRUE ~ "strong_expected_acute_innate_signal"
  )
}

classify_alignment <- function(dataset_id, manuscript_group, delta, time_h, scoring_status) {
  dataset_id <- as_chr(dataset_id)
  manuscript_group <- as_chr(manuscript_group)
  scoring_status <- as_chr(scoring_status)
  case_when(
    scoring_status != "scored_in_current_transfer_eval" & str_detect(dataset_id, "GSE237068|GSE190850") ~ "not_scored",
    scoring_status != "scored_in_current_transfer_eval" ~ "not_scored",
    str_detect(dataset_id, "GSE237068") ~ "not_directly_comparable",
    str_detect(dataset_id, "GSE166655|GSE262515|GSE314070") ~ "context_explains_weaker_signal",
    str_detect(dataset_id, "GSE264344") & !is.na(time_h) & time_h == 72 ~ "context_explains_weaker_signal",
    manuscript_group == "Primary acute validation" & !is.na(delta) & delta > 0 ~ "matches_paper_context",
    str_detect(dataset_id, "GSE264344") & !is.na(time_h) & time_h <= 24 & !is.na(delta) & delta > 0 ~ "matches_paper_context",
    !is.na(delta) & delta > 0 ~ "matches_paper_context",
    !is.na(delta) & delta <= 0 ~ "requires_manual_review",
    TRUE ~ "requires_manual_review"
  )
}

make_verification_status <- function(missing_publication, missing_split, missing_tc, missing_tt, conflict = FALSE) {
  mapply(
    function(pub, split, tc, tt, conflict_one) {
      issues <- c(
        if (isTRUE(conflict_one)) "metadata_conflict_requires_review",
        if (isTRUE(pub)) "needs_publication_verification",
        if (isTRUE(split)) "needs_split_metadata_verification",
        if (isTRUE(tc)) "needs_treatment_control_verification",
        if (isTRUE(tt)) "needs_tissue_time_verification"
      )
      if (length(issues) == 0) "verified_from_local_metadata_and_manual_paper_context" else paste(issues, collapse = "; ")
    },
    missing_publication, missing_split, missing_tc, missing_tt, conflict,
    USE.NAMES = FALSE
  )
}

scored_supp <- role %>%
  mutate(pub_key = pub_key_from_dataset(dataset_id, gse_id)) %>%
  left_join(eval_keep, by = "split_id", suffix = c("", "_eval")) %>%
  left_join(split_meta, by = "split_id") %>%
  left_join(weak_keep, by = "split_id", suffix = c("", "_weak")) %>%
  left_join(weak_interp_keep, by = "split_id", suffix = c("", "_weak_interp")) %>%
  left_join(audit_keep, by = "split_id", suffix = c("", "_audit")) %>%
  left_join(manual_paper_context_mapping, by = "pub_key") %>%
  rowwise() %>%
  mutate(source_design_file = find_design_file(dataset_id, gse_id)) %>%
  ungroup() %>%
  mutate(
    gse_id_out = coalesce_chr(gse_id_weak, gse_id, gse_id_clean, accession_manual, pub_key),
    accession_out = coalesce_chr(accession_manual, base_gse(gse_id_out), gse_id_out),
    source_database_out = coalesce_chr(source_database_manual, if_else(str_detect(coalesce_chr(accession_out, dataset_id), "^GSE"), "GEO", NA_character_)),
    PMID_out = coalesce_chr(PMID_clean, PMID, PMID_manual),
    DOI_out = coalesce_chr(DOI_clean, DOI, DOI_manual),
    paper_title_out = coalesce_chr(paper_title_clean, paper_title, paper_title_manual),
    publication_year_out = suppressWarnings(as.integer(coalesce(publication_year, publication_year_manual))),
    organism_out = coalesce_chr(organism_clean, organism_meta, organism, organism_audit, organism_original, organism_manual),
    tissue_out = coalesce_chr(tissue_meta, tissue, tissue_weak, parse_tissue(split_id, source_split_file)),
    time_h_out = suppressWarnings(as.numeric(coalesce(time_h_meta, time_h, time_h_weak, parse_time(split_id, source_split_file)))),
    timepoint_original_out = coalesce_chr(timepoint_original, ifelse(!is.na(time_h_out), paste0(time_h_out, " h"), NA_character_)),
    cell_type_or_cell_line_out = coalesce_chr(cell_line_meta, if_else(str_detect(str_to_lower(coalesce_chr(tissue_out, dataset_id)), "cell.line|cell line|hepg2"), "cell line", NA_character_)),
    in_vivo_or_in_vitro_out = case_when(
      str_detect(str_to_lower(coalesce_chr(tissue_out, dataset_id, cell_type_or_cell_line_out)), "cell.line|cell line|hepg2") ~ "in vitro",
      !is.na(tissue_out) ~ "in vivo",
      TRUE ~ NA_character_
    ),
    disease_model_out = case_when(
      str_detect(dataset_id, "GSE166655") ~ "SBMA therapeutic transgene/disease-rescue model",
      str_detect(dataset_id, "GSE262515") ~ "HBV therapeutic RNAi context",
      str_detect(dataset_id, "GSE314070") ~ "rotavirus vaccine mucosal immunity context",
      str_detect(str_to_lower(coalesce_chr(general_biology, explanation_category)), "disease|rescue|infection|rotavirus|hepatitis") ~ coalesce_chr(general_biology, explanation_category),
      TRUE ~ NA_character_
    ),
    treatment_raw_from_split = parse_token(split_id, "G"),
    control_raw_from_split = parse_token(split_id, "VS"),
    treatment_group_raw_out = coalesce_chr(treatment_group_raw_meta, treatment_group, treatment_raw_from_split, contrast_label, delivery_platform_original),
    control_group_raw_out = coalesce_chr(control_group_raw_meta, control_group, control_raw_from_split, control_definition, control_label),
    delivery_platform_out = coalesce_chr(delivery_platform, platform_meta, delivery_platform_clean, delivery_platform_type, delivery_platform_original, delivery_platform_manual),
    vector_or_formulation_out = coalesce_chr(vector_meta, formulation_meta, delivery_platform_manual, delivery_platform, delivery_platform_original, delivery_platform_clean),
    cargo_or_transgene_out = coalesce_chr(cargo_or_transgene, cargo_meta, transgene_meta),
    route_out = coalesce_chr(route_meta),
    dose_out = coalesce_chr(dose_meta),
    pipeline_role_out = case_when(
      production_anchor_row %in% TRUE ~ "Anchor",
      production_discovery_row %in% TRUE ~ "Anchor",
      original_dataset_type == "anchor" ~ "anchor_phase",
      TRUE ~ "validation"
    ),
    manuscript_group_out = coalesce_chr(final_display_group_v2,
                                        unname(role_labels[as.character(manuscript_role)])),
    manuscript_group_out = if_else(is.na(manuscript_group_out), "Excluded/unclear", manuscript_group_out),
    used_for_primary_claim_out = manuscript_group_out == "Primary acute validation",
    n_delivery_out = suppressWarnings(as.integer(coalesce(n_delivery_meta, n_delivery, n_delivery_eval))),
    n_control_out = suppressWarnings(as.integer(coalesce(n_control_meta, n_controls, n_controls_eval))),
    n_total_out = n_delivery_out + n_control_out,
    sample_ids_delivery_out = coalesce_chr(sample_ids_delivery, delivery_sample_ids),
    sample_ids_control_out = coalesce_chr(sample_ids_control, control_sample_ids),
    condition_simple_out = condition_simple_meta,
    delta_out = suppressWarnings(as.numeric(coalesce(delta_mean_imrs_z, delta_mean_imrs_z_eval, original_IMRS_delta))),
    auc_out = suppressWarnings(as.numeric(coalesce(auc_imrs_z, auc_imrs_z_eval, original_AUC_secondary))),
    direction_positive_out = if_else(is.na(delta_out), NA, delta_out > 0),
    qc_pass_out = case_when(
      pass %in% TRUE ~ TRUE,
      pass %in% FALSE ~ FALSE,
      as.character(pass) == "TRUE" ~ TRUE,
      as.character(pass) == "FALSE" ~ FALSE,
      TRUE ~ NA
    ),
    weak_score_explanation_out = coalesce_chr(manuscript_ready_interpretation, manuscript_ready_interpretation_weak_interp, explanation_category),
    reviewer_risk_out = coalesce_chr(reviewer_risk, reviewer_risk_level, reviewer_risk_weak_interp, reviewer_risk_level_weak_interp),
    notes_out = coalesce_chr(notes, notes_weak_interp, role_note, "Fields left as NA were not verified in local files. AUC is secondary."),
    scoring_status = "scored_in_current_transfer_eval",
    imrs_result_interpretation = classify_interpretation(dataset_id, manuscript_group_out, delta_out, time_h_out, expected_imrs_behavior),
    paper_result_alignment = classify_alignment(dataset_id, manuscript_group_out, delta_out, time_h_out, scoring_status),
    publication_missing = is.na(PMID_out) | is.na(paper_title_out),
    split_missing = is.na(source_split_file) | !file.exists(source_split_file),
    treatment_control_missing = is.na(treatment_group_raw_out) | is.na(control_group_raw_out),
    tissue_time_missing = is.na(tissue_out) | is.na(time_h_out),
    local_manual_organism_conflict = !is.na(organism_meta) & !is.na(organism_manual) & !str_detect(organism_manual, fixed(organism_meta)),
    verification_status = make_verification_status(
      publication_missing, split_missing, treatment_control_missing, tissue_time_missing,
      local_manual_organism_conflict
    ),
    paper_context_source = "manual_curated_from_paper"
  ) %>%
  transmute(
    dataset_id,
    gse_id = gse_id_out,
    split_id,
    source_design_file,
    source_split_file,
    sample_ids_delivery = sample_ids_delivery_out,
    sample_ids_control = sample_ids_control_out,
    n_delivery = n_delivery_out,
    n_control = n_control_out,
    n_total = n_total_out,
    organism = organism_out,
    tissue = tissue_out,
    cell_type_or_cell_line = cell_type_or_cell_line_out,
    in_vivo_or_in_vitro = in_vivo_or_in_vitro_out,
    disease_model = disease_model_out,
    time_h = time_h_out,
    timepoint_original = timepoint_original_out,
    treatment_group_raw = treatment_group_raw_out,
    control_group_raw = control_group_raw_out,
    treatment_group_clean = clean_label(treatment_group_raw_out),
    control_group_clean = clean_label(control_group_raw_out),
    condition_simple = condition_simple_out,
    delivery_platform = delivery_platform_out,
    vector_or_formulation = vector_or_formulation_out,
    cargo_or_transgene = cargo_or_transgene_out,
    route = route_out,
    dose = dose_out,
    PMID = PMID_out,
    DOI = DOI_out,
    paper_title = paper_title_out,
    publication_year = publication_year_out,
    source_database = source_database_out,
    accession = accession_out,
    pipeline_role = pipeline_role_out,
    manuscript_interpretation_group = manuscript_group_out,
    used_for_primary_claim = used_for_primary_claim_out,
    delta_mean_imrs_z = delta_out,
    auc_imrs_z_secondary = auc_out,
    direction_positive = direction_positive_out,
    qc_pass = qc_pass_out,
    inclusion_rationale,
    limitation,
    weak_score_explanation = weak_score_explanation_out,
    reviewer_risk = reviewer_risk_out,
    notes = notes_out,
    verification_status,
    scoring_status,
    expected_imrs_behavior,
    paper_result_alignment_rule,
    imrs_result_interpretation,
    paper_result_alignment,
    manuscript_context_note,
    limitations_note,
    paper_treatment_context,
    paper_context_source
  )

scored_base_keys_for_split_filter <- unique(pub_key_from_dataset(scored_supp$dataset_id, scored_supp$gse_id))
unscored_split_ids <- setdiff(split_inventory$split_id, scored_supp$split_id)
metadata_split_rows <- split_inventory %>%
  filter(split_id %in% unscored_split_ids) %>%
  filter(!(dataset_id_guess == "GSE262515" & "GSE262515" %in% scored_base_keys_for_split_filter)) %>%
  mutate(
    dataset_id = dataset_id_guess,
    pub_key = pub_key_from_dataset(dataset_id, gse_id)
  ) %>%
  left_join(manual_paper_context_mapping, by = "pub_key") %>%
  rowwise() %>%
  mutate(source_design_file = find_design_file(dataset_id, gse_id)) %>%
  ungroup() %>%
  transmute(
    dataset_id,
    gse_id = coalesce_chr(gse_id, accession_manual, pub_key),
    split_id,
    source_design_file,
    source_split_file,
    sample_ids_delivery = NA_character_,
    sample_ids_control = NA_character_,
    n_delivery,
    n_control,
    n_total = n_delivery + n_control,
    organism = organism_manual,
    tissue = coalesce_chr(tissue_recovered, tissue_from_filename),
    cell_type_or_cell_line = if_else(str_detect(str_to_lower(coalesce_chr(tissue_recovered, dataset_id)), "cell.line|cell line|balf"), coalesce_chr(tissue_recovered), NA_character_),
    in_vivo_or_in_vitro = case_when(
      str_detect(str_to_lower(coalesce_chr(tissue_recovered, dataset_id)), "cell.line|cell line") ~ "in vitro",
      TRUE ~ "in vivo"
    ),
    disease_model = NA_character_,
    time_h = coalesce(time_h_recovered, time_h_from_filename),
    timepoint_original = if_else(!is.na(time_h), paste0(time_h, " h"), NA_character_),
    treatment_group_raw,
    control_group_raw,
    treatment_group_clean = clean_label(treatment_group_raw),
    control_group_clean = clean_label(control_group_raw),
    condition_simple = NA_character_,
    delivery_platform = delivery_platform_manual,
    vector_or_formulation = delivery_platform_manual,
    cargo_or_transgene = NA_character_,
    route = NA_character_,
    dose = NA_character_,
    PMID = PMID_manual,
    DOI = DOI_manual,
    paper_title = paper_title_manual,
    publication_year = publication_year_manual,
    source_database = source_database_manual,
    accession = accession_manual,
    pipeline_role = "validation",
    manuscript_interpretation_group = "Excluded/unclear",
    used_for_primary_claim = FALSE,
    delta_mean_imrs_z = NA_real_,
    auc_imrs_z_secondary = NA_real_,
    direction_positive = NA,
    qc_pass = NA,
    inclusion_rationale = "Verified metadata and split design exist, but this split is not scored in the current IMRS transfer evaluation.",
    limitation = limitations_note,
    weak_score_explanation = NA_character_,
    reviewer_risk = "manual review if later used",
    notes = "Metadata-only split. No current IMRS result is assigned.",
    verification_status = make_verification_status(is.na(PMID_manual) | is.na(paper_title_manual), FALSE, is.na(treatment_group_raw) | is.na(control_group_raw), is.na(tissue) | is.na(time_h)),
    scoring_status = "verified_metadata_only_not_scored",
    expected_imrs_behavior,
    paper_result_alignment_rule,
    imrs_result_interpretation = classify_interpretation(dataset_id, manuscript_interpretation_group, delta_mean_imrs_z, time_h, expected_imrs_behavior),
    paper_result_alignment = classify_alignment(dataset_id, manuscript_interpretation_group, delta_mean_imrs_z, time_h, scoring_status),
    manuscript_context_note,
    limitations_note,
    paper_treatment_context,
    paper_context_source = "manual_curated_from_paper"
  )

scored_base_keys <- unique(pub_key_from_dataset(scored_supp$dataset_id, scored_supp$gse_id))
metadata_split_base_keys <- unique(pub_key_from_dataset(metadata_split_rows$dataset_id, metadata_split_rows$gse_id))
design_only_expected <- setdiff(
  expected_metadata_folders,
  c(scored_supp$dataset_id, metadata_split_rows$dataset_id, scored_base_keys, metadata_split_base_keys, "GSE262515")
)

design_only_rows <- tibble(dataset_id = design_only_expected) %>%
  filter(!is.na(dataset_id), dataset_id != "") %>%
  mutate(
    gse_id = base_gse(dataset_id),
    pub_key = pub_key_from_dataset(dataset_id, gse_id)
  ) %>%
  left_join(manual_paper_context_mapping, by = "pub_key") %>%
  rowwise() %>%
  mutate(source_design_file = find_design_file(dataset_id, gse_id)) %>%
  ungroup() %>%
  transmute(
    dataset_id,
    gse_id = coalesce_chr(dataset_id, accession_manual),
    split_id = paste0(dataset_id, "_metadata_only"),
    source_design_file,
    source_split_file = NA_character_,
    sample_ids_delivery = NA_character_,
    sample_ids_control = NA_character_,
    n_delivery = NA_integer_,
    n_control = NA_integer_,
    n_total = NA_integer_,
    organism = organism_manual,
    tissue = if_else(dataset_id == "GSE237068", "Liver", NA_character_),
    cell_type_or_cell_line = NA_character_,
    in_vivo_or_in_vitro = if_else(dataset_id == "GSE237068", "in vivo", NA_character_),
    disease_model = NA_character_,
    time_h = if_else(dataset_id == "GSE237068", 96, NA_real_),
    timepoint_original = if_else(dataset_id == "GSE237068", "4 days after retro-orbital adenovirus administration", NA_character_),
    treatment_group_raw = if_else(dataset_id == "GSE237068", "delivery_ecsth", NA_character_),
    control_group_raw = if_else(dataset_id == "GSE237068", "delivery_gfp", NA_character_),
    treatment_group_clean = clean_label(treatment_group_raw),
    control_group_clean = clean_label(control_group_raw),
    condition_simple = if_else(dataset_id == "GSE237068", "DELIVERY", NA_character_),
    delivery_platform = delivery_platform_manual,
    vector_or_formulation = if_else(dataset_id == "GSE237068", "EcSTH-expressing adenovirus versus GFP adenovirus", delivery_platform_manual),
    cargo_or_transgene = if_else(dataset_id == "GSE237068", "EcSTH; GFP control", NA_character_),
    route = if_else(dataset_id == "GSE237068", "retro-orbital", NA_character_),
    dose = NA_character_,
    PMID = PMID_manual,
    DOI = DOI_manual,
    paper_title = paper_title_manual,
    publication_year = publication_year_manual,
    source_database = source_database_manual,
    accession = accession_manual,
    pipeline_role = "validation",
    manuscript_interpretation_group = "Excluded/unclear",
    used_for_primary_claim = FALSE,
    delta_mean_imrs_z = NA_real_,
    auc_imrs_z_secondary = NA_real_,
    direction_positive = NA,
    qc_pass = NA,
    inclusion_rationale = if_else(dataset_id == "GSE237068", "Verified liver adenoviral-vector metadata retained for provenance; excluded from current scoring because the EcSTH metabolic perturbation is not a primary delivery-response experiment.", "Verified metadata folder/design exists, but no current scored split is present in the manuscript role table."),
    limitation = limitations_note,
    weak_score_explanation = NA_character_,
    reviewer_risk = "manual review if later used",
    notes = if_else(dataset_id == "GSE237068", "Metadata-only dataset; liver collected 96 h after retro-orbital EcSTH- or GFP-adenovirus administration. No current IMRS result is assigned.", "Metadata-only dataset. No current IMRS result is assigned."),
    verification_status = if_else(dataset_id == "GSE237068", "verified_metadata_only_not_scored; treatment_control_tissue_time_verified", make_verification_status(is.na(PMID_manual) | is.na(paper_title_manual), is.na(source_design_file), TRUE, TRUE)),
    scoring_status = "verified_metadata_only_not_scored",
    expected_imrs_behavior,
    paper_result_alignment_rule,
    imrs_result_interpretation = classify_interpretation(dataset_id, manuscript_interpretation_group, delta_mean_imrs_z, time_h, expected_imrs_behavior),
    paper_result_alignment = classify_alignment(dataset_id, manuscript_interpretation_group, delta_mean_imrs_z, time_h, scoring_status),
    manuscript_context_note,
    limitations_note,
    paper_treatment_context,
    paper_context_source = "manual_curated_from_paper"
  )

supp <- bind_rows(scored_supp, metadata_split_rows, design_only_rows) %>%
  mutate(
    manuscript_interpretation_group = factor(
      manuscript_interpretation_group,
      levels = c("Locked anchor", "Primary acute validation", "Extended validation",
                 "Secondary support", "Excluded/unclear")
    ),
    manuscript_interpretation_group = as.character(manuscript_interpretation_group),
    used_for_primary_claim = manuscript_interpretation_group == "Primary acute validation",
    notes = coalesce_chr(notes, "Fields left as NA were not verified in local files. AUC is secondary.")
  ) %>%
  arrange(scoring_status != "scored_in_current_transfer_eval", dataset_id, time_h, tissue, split_id)

folder_inventory <- tibble(expected_dataset_folder = expected_metadata_folders) %>%
  mutate(
    base_key = pub_key_from_dataset(expected_dataset_folder, expected_dataset_folder),
    detected_top_level_design_file = vapply(expected_dataset_folder, find_design_file, character(1)),
    detected_scoring_design_files = vapply(expected_dataset_folder, function(x) {
      key <- base_gse(x)
      hits <- design_files %>% filter(under_scoring_dir, design_key == x | dataset_key_base == key)
      collapse_unique(hits$path_norm, max_items = 5)
    }, character(1)),
    detected_split_directories = vapply(expected_dataset_folder, function(x) {
      key <- base_gse(x)
      hits <- split_file_index %>% filter(dataset_id_guess == x | gse_id == key | base_gse(dataset_id_guess) == key)
      collapse_unique(dirname(hits$source_split_file), max_items = 5)
    }, character(1)),
    n_split_design_files = vapply(expected_dataset_folder, function(x) {
      key <- base_gse(x)
      nrow(split_file_index %>% filter(dataset_id_guess == x | gse_id == key | base_gse(dataset_id_guess) == key))
    }, integer(1)),
    n_scored_rows = vapply(expected_dataset_folder, function(x) {
      key <- base_gse(x)
      nrow(scored_supp %>% filter(dataset_id == x | base_gse(dataset_id) == key))
    }, integer(1)),
    has_verified_metadata = !is.na(detected_top_level_design_file) | n_split_design_files > 0 | !is.na(detected_scoring_design_files),
    scoring_status_summary = case_when(
      n_scored_rows > 0 ~ "scored_in_current_transfer_eval",
      has_verified_metadata ~ "verified_metadata_only_not_scored",
      TRUE ~ "missing_from_verified_metadata"
    ),
    notes = case_when(
      expected_dataset_folder == "GSE262515" ~ "Scored rows are represented by GSE262515_cell_line and GSE262515_tissue.",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-base_key)

write_tsv(folder_inventory, out_paths$folder_inventory, na = "NA")

interpretation_group_summary <- supp %>%
  group_by(scoring_status, manuscript_interpretation_group) %>%
  summarise(
    split_rows = n(),
    datasets = n_distinct(dataset_id),
    primary_claim_rows = sum(used_for_primary_claim %in% TRUE, na.rm = TRUE),
    mean_delta_mean_imrs_z = if_else(all(is.na(delta_mean_imrs_z)), NA_real_, mean(delta_mean_imrs_z, na.rm = TRUE)),
    positive_directionality = if_else(all(is.na(direction_positive)), NA_real_, mean(direction_positive %in% TRUE, na.rm = TRUE)),
    .groups = "drop"
  )

dataset_summary <- supp %>%
  group_by(dataset_id, scoring_status, manuscript_interpretation_group, PMID, paper_title) %>%
  summarise(
    split_rows = n(),
    tissues = collapse_unique(tissue, max_items = 5),
    time_h_range = if_else(all(is.na(time_h)), NA_character_, paste0(min(time_h, na.rm = TRUE), " to ", max(time_h, na.rm = TRUE))),
    delivery_platforms = collapse_unique(delivery_platform, max_items = 4),
    mean_delta_mean_imrs_z = if_else(all(is.na(delta_mean_imrs_z)), NA_real_, mean(delta_mean_imrs_z, na.rm = TRUE)),
    paper_result_alignment = collapse_unique(paper_result_alignment, max_items = 4),
    .groups = "drop"
  )

publication_summary <- supp %>%
  group_by(PMID, DOI, paper_title, publication_year) %>%
  summarise(
    datasets = collapse_unique(dataset_id, max_items = 6),
    rows = n(),
    scoring_status = collapse_unique(scoring_status, max_items = 4),
    manuscript_groups = collapse_unique(manuscript_interpretation_group, max_items = 5),
    .groups = "drop"
  ) %>%
  arrange(is.na(PMID), PMID)

treatment_control_summary <- supp %>%
  group_by(dataset_id, treatment_group_clean, control_group_clean, delivery_platform, tissue) %>%
  summarise(
    rows = n(),
    time_h_range = if_else(all(is.na(time_h)), NA_character_, paste0(min(time_h, na.rm = TRUE), " to ", max(time_h, na.rm = TRUE))),
    mean_delta_mean_imrs_z = if_else(all(is.na(delta_mean_imrs_z)), NA_real_, mean(delta_mean_imrs_z, na.rm = TRUE)),
    paper_result_alignment = collapse_unique(paper_result_alignment, max_items = 3),
    .groups = "drop"
  ) %>%
  arrange(dataset_id, tissue, treatment_group_clean)

imrs_result_alignment_summary <- supp %>%
  group_by(scoring_status, imrs_result_interpretation, paper_result_alignment) %>%
  summarise(rows = n(), datasets = n_distinct(dataset_id), .groups = "drop") %>%
  arrange(scoring_status, paper_result_alignment)

missing_or_needs_verification <- supp %>%
  filter(
    verification_status != "verified_from_local_metadata_and_manual_paper_context" |
      is.na(PMID) | is.na(paper_title) |
      is.na(treatment_group_raw) | is.na(control_group_raw) |
      is.na(tissue) | is.na(time_h)
  ) %>%
  transmute(
    dataset_id,
    split_id,
    scoring_status,
    manuscript_interpretation_group,
    verification_status,
    missing_PMID = is.na(PMID),
    missing_title = is.na(paper_title),
    missing_treatment_control = is.na(treatment_group_raw) | is.na(control_group_raw),
    missing_tissue_time = is.na(tissue) | is.na(time_h),
    paper_result_alignment
  )

weak_dataset_context <- supp %>%
  filter(
    paper_result_alignment == "context_explains_weaker_signal" |
      imrs_result_interpretation == "weak_but_biologically_explainable" |
      str_detect(dataset_id, "GSE166655|GSE262515|GSE264344|GSE314070")
  ) %>%
  select(
    dataset_id, split_id, tissue, time_h, treatment_group_clean, control_group_clean,
    delta_mean_imrs_z, manuscript_interpretation_group, expected_imrs_behavior,
    imrs_result_interpretation, paper_result_alignment, paper_treatment_context,
    paper_result_alignment_rule, limitations_note
  )

readme_sheet <- tibble(
  field = c(
    "description",
    "row_definition",
    "scored_rows",
    "metadata_only_policy",
    "publication_policy",
    "primary_claim_policy",
    "auc_policy"
  ),
  value = c(
    "All-dataset, all-split provenance and biological-context audit for the IMRS project.",
    "One row per scored split contrast, plus rows for verified metadata/split designs that are not scored in the current transfer evaluation.",
    "Scored rows are sourced from manuscript_dataset_role_table.tsv and split-level IMRS evaluation outputs.",
    "Metadata-only rows retain verified metadata and paper context but do not assign IMRS results.",
    "Publication fields are populated from manual curated paper-context mapping and local audit files; unverified fields are NA.",
    "used_for_primary_claim is TRUE only for Primary acute validation.",
    "AUC is retained only as a secondary metric."
  )
)

write_tsv(supp, out_paths$tsv, na = "NA")
write_csv(supp, out_paths$csv, na = "NA")

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      README = readme_sheet,
      folder_inventory_all_verified = folder_inventory,
      split_design_inventory_all = split_inventory,
      full_split_provenance_all = supp,
      scored_split_provenance_only = supp %>% filter(scoring_status == "scored_in_current_transfer_eval"),
      metadata_only_not_scored = supp %>% filter(scoring_status != "scored_in_current_transfer_eval"),
      dataset_summary = dataset_summary,
      interpretation_group_summary = interpretation_group_summary,
      publication_summary = publication_summary,
      treatment_control_summary = treatment_control_summary,
      imrs_result_alignment_summary = imrs_result_alignment_summary,
      missing_or_needs_verification = missing_or_needs_verification,
      weak_dataset_context = weak_dataset_context,
      manual_paper_context_mapping = manual_paper_context_mapping
    ),
    out_paths$xlsx
  )
} else {
  warning("writexl is not installed; XLSX output was not written.")
}

supp_rmd <- c(
  "---",
  "title: \"Supplementary Dataset and Split-Level Provenance Audit\"",
  "subtitle: \"All-dataset biological-context audit for the IMRS publication-readiness assessment\"",
  "author: \"Yansheng Luo\"",
  "date: \"`r format(Sys.Date(), '%B %d, %Y')`\"",
  "output:",
  "  pdf_document:",
  "    toc: true",
  "    toc_depth: 2",
  "    number_sections: true",
  "    latex_engine: xelatex",
  "fontsize: 10pt",
  "geometry: margin=0.72in",
  "header-includes:",
  "  - \\usepackage{booktabs}",
  "  - \\usepackage{array}",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "library(readr); library(dplyr); library(knitr); library(stringr)",
  paste0("supp <- readr::read_tsv('", out_paths$tsv, "', show_col_types = FALSE, na = c('', 'NA'))"),
  paste0("folder_inventory <- readr::read_tsv('", out_paths$folder_inventory, "', show_col_types = FALSE, na = c('', 'NA'))"),
  paste0("split_inventory <- readr::read_tsv('", out_paths$split_inventory, "', show_col_types = FALSE, na = c('', 'NA'))"),
  "fmt_num <- function(x, digits = 3) ifelse(is.na(x), 'NA', formatC(as.numeric(x), digits = digits, format = 'f'))",
  "fmt_pct <- function(x) ifelse(is.na(x), 'NA', paste0(formatC(as.numeric(x) * 100, digits = 1, format = 'f'), '%'))",
  "short <- function(x, n = 70) stringr::str_trunc(as.character(x), n)",
  "kable_compact <- function(x, caption) knitr::kable(x, format = 'latex', booktabs = TRUE, longtable = FALSE, caption = caption, escape = TRUE)",
  "```",
  "",
  "# Purpose",
  "",
  "This supplement provides all-dataset and all-split provenance for the IMRS publication-readiness assessment. It links verified metadata, split-design files, current IMRS scoring results, source publication context, treatment/control design, tissue, timepoint, delivery platform, and expected biological interpretation.",
  "",
  "The complete table is provided as machine-readable XLSX/TSV/CSV because the full split-level provenance table is too wide for a readable PDF page.",
  "",
  "# Data Sources Used",
  "",
  "```{r sources}",
  "sources <- tibble::tibble(",
  "  Source = c('Cleaned manuscript role table', 'Split-level IMRS evaluation', 'Verified metadata directory', 'Weak/context audit tables', 'Manual curated paper-context mapping'),",
  "  Use = c('Defines scored split contrasts and manuscript interpretation groups', 'Adds Delta IMRSz, AUC as secondary metric, QC and sample counts', 'Recovers source design files, split files, sample IDs, treatment/control labels, tissue and timepoint', 'Adds weak-score interpretation and reviewer-risk context where available', 'Adds PMID, DOI, paper title, expected IMRS behavior and paper-result alignment rules')",
  ")",
  "kable_compact(sources, 'Data sources used to create the supplement.')",
  "```",
  "",
  "# All Verified Metadata Folders Audited",
  "",
  "```{r folder-inventory}",
  "folder_preview <- folder_inventory %>%",
  "  transmute(",
  "    dataset = expected_dataset_folder,",
  "    split_files = n_split_design_files,",
  "    scored_rows = n_scored_rows,",
  "    scoring_status = scoring_status_summary,",
  "    metadata_detected = has_verified_metadata,",
  "    note = short(notes, 42)",
  "  )",
  "kable_compact(folder_preview, 'Verified metadata folders audited. Full file paths are stored in the XLSX and TSV inventories.')",
  "```",
  "",
  "# Summary by Scoring Status",
  "",
  "```{r scoring-status}",
  "scoring_summary <- supp %>%",
  "  group_by(scoring_status) %>%",
  "  summarise(rows = n(), datasets = n_distinct(dataset_id), publication_metadata_rows = sum(!is.na(PMID) & !is.na(paper_title)), .groups = 'drop')",
  "kable_compact(scoring_summary, 'Summary by scoring status.')",
  "```",
  "",
  "# Summary by Manuscript Interpretation Group",
  "",
  "```{r group-summary}",
  "group_summary <- supp %>%",
  "  group_by(manuscript_interpretation_group) %>%",
  "  summarise(",
  "    rows = n(),",
  "    datasets = n_distinct(dataset_id),",
  "    primary_claim_rows = sum(used_for_primary_claim %in% TRUE, na.rm = TRUE),",
  "    mean_delta = fmt_num(ifelse(all(is.na(delta_mean_imrs_z)), NA, mean(delta_mean_imrs_z, na.rm = TRUE)), 3),",
  "    positive_directionality = fmt_pct(ifelse(all(is.na(direction_positive)), NA, mean(direction_positive %in% TRUE, na.rm = TRUE))),",
  "    .groups = 'drop'",
  "  )",
  "kable_compact(group_summary, 'Summary by manuscript interpretation group.')",
  "```",
  "",
  "# Summary by Dataset",
  "",
  "```{r dataset-summary}",
  "dataset_summary <- supp %>%",
  "  group_by(dataset_id, scoring_status, manuscript_interpretation_group) %>%",
  "  summarise(",
  "    rows = n(),",
  "    tissues = short(paste(sort(unique(na.omit(tissue))), collapse = '; '), 38),",
  "    time_h_min = suppressWarnings(min(time_h, na.rm = TRUE)),",
  "    time_h_max = suppressWarnings(max(time_h, na.rm = TRUE)),",
  "    mean_delta = fmt_num(ifelse(all(is.na(delta_mean_imrs_z)), NA, mean(delta_mean_imrs_z, na.rm = TRUE)), 3),",
  "    alignment = short(paste(sort(unique(na.omit(paper_result_alignment))), collapse = '; '), 42),",
  "    .groups = 'drop'",
  "  ) %>%",
  "  mutate(time_h_min = ifelse(is.infinite(time_h_min), NA, time_h_min), time_h_max = ifelse(is.infinite(time_h_max), NA, time_h_max))",
  "kable_compact(dataset_summary, 'Summary by dataset.')",
  "```",
  "",
  "# Summary by Publication Source",
  "",
  "```{r publication-summary}",
  "publication_summary <- supp %>%",
  "  group_by(PMID, DOI, paper_title) %>%",
  "  summarise(",
  "    datasets = short(paste(sort(unique(dataset_id)), collapse = '; '), 44),",
  "    rows = n(),",
  "    scoring = short(paste(sort(unique(scoring_status)), collapse = '; '), 34),",
  "    .groups = 'drop'",
  "  ) %>%",
  "  mutate(paper_title = short(paper_title, 58)) %>%",
  "  arrange(is.na(PMID), PMID)",
  "kable_compact(publication_summary, 'Summary by source publication / PMID.')",
  "```",
  "",
  "# Summary by Treatment/Control Design",
  "",
  "```{r treatment-summary}",
  "tc_summary <- supp %>%",
  "  filter(scoring_status == 'scored_in_current_transfer_eval') %>%",
  "  group_by(dataset_id, tissue, treatment_group_clean, control_group_clean) %>%",
  "  summarise(rows = n(), time_h_range = paste0(suppressWarnings(min(time_h, na.rm = TRUE)), ' to ', suppressWarnings(max(time_h, na.rm = TRUE))), mean_delta = fmt_num(mean(delta_mean_imrs_z, na.rm = TRUE), 3), .groups = 'drop') %>%",
  "  mutate(treatment_group_clean = short(treatment_group_clean, 28), control_group_clean = short(control_group_clean, 24), tissue = short(tissue, 20), time_h_range = ifelse(str_detect(time_h_range, 'Inf'), 'NA', time_h_range)) %>%",
  "  slice_head(n = 18)",
  "kable_compact(tc_summary, 'Compact preview of scored treatment/control designs. Complete labels are in the XLSX/TSV/CSV.')",
  "```",
  "",
  "# IMRS Result Alignment with Paper Context",
  "",
  "```{r alignment-summary}",
  "alignment_summary <- supp %>%",
  "  group_by(imrs_result_interpretation, paper_result_alignment) %>%",
  "  summarise(rows = n(), datasets = n_distinct(dataset_id), .groups = 'drop')",
  "kable_compact(alignment_summary, 'IMRS result interpretation and paper-context alignment summary.')",
  "```",
  "",
  "```{r compact-context-table}",
  "compact_context <- supp %>%",
  "  group_by(dataset_id, scoring_status, manuscript_interpretation_group, paper_treatment_context, imrs_result_interpretation, paper_result_alignment) %>%",
  "  summarise(rows = n(), mean_delta = ifelse(all(is.na(delta_mean_imrs_z)), NA, mean(delta_mean_imrs_z, na.rm = TRUE)), .groups = 'drop') %>%",
  "  mutate(",
  "    scored = ifelse(scoring_status == 'scored_in_current_transfer_eval', 'yes', 'no'),",
  "    `paper context` = short(paper_treatment_context, 56),",
  "    `observed IMRS behavior` = ifelse(is.na(mean_delta), imrs_result_interpretation, paste0(imrs_result_interpretation, '; mean Delta=', fmt_num(mean_delta, 2))),",
  "    `paper-result alignment` = paper_result_alignment",
  "  ) %>%",
  "  select(dataset_id, scored, `manuscript group` = manuscript_interpretation_group, `paper context`, `observed IMRS behavior`, `paper-result alignment`) %>%",
  "  arrange(dataset_id)",
  "kable_compact(compact_context, 'Compact dataset-level paper-context interpretation table. Full wording and split-level details are in the XLSX/TSV/CSV.')",
  "```",
  "",
  "# Metadata-Only Datasets Not Currently Scored",
  "",
  "```{r metadata-only}",
  "metadata_only <- supp %>%",
  "  filter(scoring_status != 'scored_in_current_transfer_eval') %>%",
  "  transmute(dataset_id, split_id = short(split_id, 42), tissue, time_h, PMID, paper_result_alignment, note = short(manuscript_context_note, 58))",
  "if (nrow(metadata_only) == 0) metadata_only <- tibble::tibble(Note = 'No metadata-only rows were detected.')",
  "kable_compact(metadata_only, 'Verified metadata rows that are not scored in the current transfer evaluation.')",
  "```",
  "",
  "# Weak-Score/Context-Shifted Dataset Provenance",
  "",
  "```{r weak-context}",
  "weak_summary <- supp %>%",
  "  filter(paper_result_alignment == 'context_explains_weaker_signal' | imrs_result_interpretation == 'weak_but_biologically_explainable') %>%",
  "  group_by(dataset_id, PMID, manuscript_interpretation_group, paper_result_alignment) %>%",
  "  summarise(rows = n(), time_h_range = paste0(suppressWarnings(min(time_h, na.rm = TRUE)), ' to ', suppressWarnings(max(time_h, na.rm = TRUE))), mean_delta = fmt_num(ifelse(all(is.na(delta_mean_imrs_z)), NA, mean(delta_mean_imrs_z, na.rm = TRUE)), 3), explanation = short(first(paper_result_alignment_rule), 66), .groups = 'drop') %>%",
  "  mutate(time_h_range = ifelse(str_detect(time_h_range, 'Inf'), 'NA', time_h_range))",
  "kable_compact(weak_summary, 'Weak-score and context-shifted dataset provenance summary.')",
  "```",
  "",
  "# Missing Fields and Verification Needs",
  "",
  "```{r missing-summary}",
  "missing_summary <- supp %>%",
  "  filter(verification_status != 'verified_from_local_metadata_and_manual_paper_context' | is.na(PMID) | is.na(paper_title)) %>%",
  "  group_by(dataset_id, verification_status, paper_result_alignment) %>%",
  "  summarise(rows = n(), missing_PMID_or_title = sum(is.na(PMID) | is.na(paper_title)), .groups = 'drop') %>%",
  "  mutate(verification_status = short(verification_status, 70))",
  "if (nrow(missing_summary) == 0) missing_summary <- tibble::tibble(Note = 'No publication-metadata gaps were detected. Metadata-only rows retain NA split fields where no scored contrast exists.')",
  "kable_compact(missing_summary, 'Rows requiring metadata or publication verification.')",
  "```",
  "",
  "# Preview of Full Provenance Table",
  "",
  "```{r preview}",
  "preview <- supp %>%",
  "  select(dataset_id, tissue, time_h, treatment_group_clean, control_group_clean, manuscript_interpretation_group, delta_mean_imrs_z, PMID, verification_status) %>%",
  "  mutate(treatment_group_clean = short(treatment_group_clean, 24), control_group_clean = short(control_group_clean, 22), delta_mean_imrs_z = fmt_num(delta_mean_imrs_z, 3), verification_status = short(verification_status, 38)) %>%",
  "  slice_head(n = 10)",
  "kable_compact(preview, 'Preview of the full provenance table. Complete split-level details are available as XLSX/TSV/CSV.')",
  "```",
  "",
  "# Output Files",
  "",
  paste0("- ", out_paths$xlsx),
  paste0("- ", out_paths$tsv),
  paste0("- ", out_paths$csv),
  paste0("- ", out_paths$folder_inventory),
  paste0("- ", out_paths$split_inventory),
  paste0("- ", out_paths$pdf)
)

writeLines(supp_rmd, out_paths$rmd, useBytes = TRUE)

qa <- list(
  verified_metadata_folders_detected = nrow(folder_inventory),
  scored_split_rows = sum(supp$scoring_status == "scored_in_current_transfer_eval"),
  metadata_only_rows = sum(supp$scoring_status == "verified_metadata_only_not_scored"),
  datasets_with_publication_metadata = n_distinct(supp$dataset_id[!is.na(supp$PMID) & !is.na(supp$paper_title)]),
  rows_with_paper_result_alignment = sum(!is.na(supp$paper_result_alignment)),
  rows_requiring_manual_review = sum(supp$paper_result_alignment == "requires_manual_review" | str_detect(supp$verification_status, "requires_review"), na.rm = TRUE),
  rows_missing_from_role_table = length(setdiff(role$split_id, supp$split_id)),
  rows_with_recovered_treatment_control = sum(!is.na(supp$treatment_group_raw) & !is.na(supp$control_group_raw))
)

cat("Folder inventory TSV:", out_paths$folder_inventory, "\n")
cat("Split design inventory TSV:", out_paths$split_inventory, "\n")
cat("Supplement TSV:", out_paths$tsv, "\n")
cat("Supplement CSV:", out_paths$csv, "\n")
cat("Supplement XLSX:", out_paths$xlsx, "\n")
cat("Supplement Rmd:", out_paths$rmd, "\n")
cat("Supplement PDF:", out_paths$pdf, "\n")
cat("Verified metadata folders detected:", qa$verified_metadata_folders_detected, "\n")
cat("Scored split rows:", qa$scored_split_rows, "\n")
cat("Metadata-only rows:", qa$metadata_only_rows, "\n")
cat("Datasets with publication metadata:", qa$datasets_with_publication_metadata, "\n")
cat("Rows with paper_result_alignment:", qa$rows_with_paper_result_alignment, "\n")
cat("Rows requiring manual review:", qa$rows_requiring_manual_review, "\n")
cat("Rows missing from manuscript role table:", qa$rows_missing_from_role_table, "\n")
cat("Rows with recovered treatment/control metadata:", qa$rows_with_recovered_treatment_control, "\n")
