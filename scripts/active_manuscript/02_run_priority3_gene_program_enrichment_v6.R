#!/usr/bin/env Rscript
# Priority 3 retained IMRS gene-program enrichment analysis for the reviewer release.
#
# This script does not refit IMRS, change retained genes, or modify scoring,
# validation, transfer-evaluation, or robustness results.

options(stringsAsFactors = FALSE)

imrs_detect_script_path <- function(expected_basename = NULL) {
  candidates <- character()

  # Rscript --file=... invocation.
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    candidates <- c(candidates, sub("^--file=", "", file_arg[[1L]]))
  }

  # source()/sys.source() invocation (search all active frames, not only frame 1).
  frame_paths <- unlist(lapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.null(path) || length(path) == 0L) character() else as.character(path[[1L]])
  }), use.names = FALSE)
  candidates <- c(candidates, frame_paths)

  # RStudio "Run" / "Source" invocation, including running selected lines.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    editor_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    candidates <- c(candidates, editor_path)
  }

  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (length(candidates) > 0L) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
    if (!is.null(expected_basename) && nzchar(expected_basename)) {
      matching <- candidates[basename(candidates) == expected_basename]
      if (length(matching) > 0L) candidates <- c(matching, candidates)
    }
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
  }

  # Last-resort support when the working directory is the script directory.
  if (!is.null(expected_basename) && nzchar(expected_basename)) {
    cwd_candidate <- file.path(getwd(), expected_basename)
    if (file.exists(cwd_candidate)) {
      return(normalizePath(cwd_candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

imrs_find_repo_root_bootstrap <- function(start = getwd()) {
  if (is.null(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
    start <- getwd()
  }
  current <- normalizePath(start, winslash = "/", mustWork = FALSE)
  if (file.exists(current) && !dir.exists(current)) current <- dirname(current)

  repeat {
    marker_config <- file.path(current, "config", "config_template.yml")
    marker_active <- file.path(current, "scripts", "active_manuscript", "lib", "active_config.R")
    if (file.exists(marker_config) && file.exists(marker_active)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  NA_character_
}

this_file <- imrs_detect_script_path("02_run_priority3_gene_program_enrichment_v6.R")
bootstrap_starts <- unique(c(
  if (!is.na(this_file)) dirname(this_file) else character(),
  getwd()
))
repo_root_bootstrap <- NA_character_
for (start in bootstrap_starts) {
  candidate_root <- imrs_find_repo_root_bootstrap(start)
  if (!is.na(candidate_root)) {
    repo_root_bootstrap <- candidate_root
    break
  }
}
if (is.na(repo_root_bootstrap)) {
  stop(
    "Could not locate the IMRS repository root. Keep the extracted repository structure intact and run this file with RStudio Source/Run or Rscript.",
    call. = FALSE
  )
}
active_config_helper <- file.path(
  repo_root_bootstrap, "scripts", "active_manuscript", "lib", "active_config.R"
)
source(active_config_helper)


config <- imrs_load_active_config(repo_root_bootstrap)
project_root <- imrs_project_root(config)
v6_root <- imrs_config_field_path(config, "manuscript_output_dir")
analysis_root <- imrs_config_field_path(config, "priority3_enrichment_dir")
tables_dir <- file.path(analysis_root, "tables")
figures_dir <- file.path(analysis_root, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(analysis_root, "Priority3_enrichment_analysis_log.txt")
log_lines <- character()
log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
  invisible(line)
}
write_log <- function() {
  writeLines(log_lines, log_path, useBytes = TRUE)
}
fail_with_log <- function(...) {
  log_msg("ERROR: ", ...)
  write_log()
  stop(paste(..., collapse = ""), call. = FALSE)
}

log_msg("Starting Priority 3 retained IMRS gene-program enrichment analysis.")
log_msg("Project root: ", project_root)
log_msg("Output folder: ", analysis_root)

core_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble", "ggplot2",
  "forcats", "scales", "purrr", "patchwork", "AnnotationDbi"
)
analysis_packages <- c("clusterProfiler", "org.Mm.eg.db", "ReactomePA", "msigdbr")
table_packages <- c("writexl", "svglite")

install_commands <- c(
  clusterProfiler = "BiocManager::install('clusterProfiler', ask = FALSE, update = FALSE)",
  org.Mm.eg.db = "BiocManager::install('org.Mm.eg.db', ask = FALSE, update = FALSE)",
  ReactomePA = "BiocManager::install('ReactomePA', ask = FALSE, update = FALSE)",
  msigdbr = "install.packages('msigdbr', repos = 'https://cloud.r-project.org')",
  writexl = "install.packages('writexl', repos = 'https://cloud.r-project.org')",
  svglite = "install.packages('svglite', repos = 'https://cloud.r-project.org')"
)

check_optional_package <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  log_msg("Optional package unavailable; it is not installed automatically by the reviewer runner: ",
          pkg, ". Install command: ", install_commands[[pkg]])
  FALSE
}

missing_core <- core_packages[!vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_core) > 0) {
  fail_with_log("Missing core package(s): ", paste(missing_core, collapse = ", "))
}
for (pkg in c(analysis_packages, table_packages)) {
  check_optional_package(pkg)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(purrr)
  library(patchwork)
  library(AnnotationDbi)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
delta_symbol <- intToUtf8(0x0394)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path, na = "")
}

strip_ensembl_version <- function(x) {
  stringr::str_replace(as.character(x), "\\.[0-9]+$", "")
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

detect_gene_column <- function(tbl, candidates) {
  hit <- candidates[candidates %in% names(tbl)]
  if (length(hit) == 0) {
    stop("Could not detect a gene ID column. Columns: ", paste(names(tbl), collapse = ", "), call. = FALSE)
  }
  hit[[1]]
}

gene_weights_path <- first_existing(c(
  imrs_config_field_path(config, "frozen_gene_weights"),
  file.path(v6_root, "tables", "gene_weights.tsv"),
  file.path(v6_root, "gene_weights.tsv")
))
if (is.na(gene_weights_path)) fail_with_log("Retained IMRS gene table not found.")
log_msg("Retained gene source: ", gene_weights_path)

weights_tbl <- readr::read_tsv(gene_weights_path, show_col_types = FALSE, progress = FALSE)
gene_col <- detect_gene_column(weights_tbl, c("gene", "gene_id", "ensembl_gene_id", "gene_id_clean", "ENSEMBL"))
query_tbl <- weights_tbl %>%
  mutate(
    query_order = dplyr::row_number(),
    original_gene_id = as.character(.data[[gene_col]]),
    ensembl_base = strip_ensembl_version(original_gene_id)
  ) %>%
  distinct(original_gene_id, .keep_all = TRUE)

n_query_original <- nrow(query_tbl)
log_msg("Retained query gene count loaded: ", n_query_original)
if (!identical(n_query_original, 287L)) {
  fail_with_log("Expected exactly 287 retained frozen IMRS genes, but loaded ", n_query_original, ".")
}

mapping_path <- first_existing(c(
  imrs_config_field_path(config, "gene_symbol_mapping"),
  file.path(v6_root, "tables", "gene_symbol_mapping.tsv"),
  file.path(v6_root, "gene_symbol_mapping.tsv")
))
project_map <- tibble(ensembl_base = character(), project_mgi_symbol = character())
if (!is.na(mapping_path)) {
  log_msg("Project gene symbol mapping source: ", mapping_path)
  raw_project_map <- readr::read_tsv(mapping_path, show_col_types = FALSE, progress = FALSE)
  if (all(c("ensembl_gene_id", "mgi_symbol") %in% names(raw_project_map))) {
    project_map <- raw_project_map %>%
      transmute(
        ensembl_base = strip_ensembl_version(.data$ensembl_gene_id),
        project_mgi_symbol = as.character(.data$mgi_symbol)
      ) %>%
      group_by(ensembl_base) %>%
      summarise(project_mgi_symbol = paste(sort(unique(na.omit(project_mgi_symbol))), collapse = ";"), .groups = "drop")
  } else {
    log_msg("Project mapping file did not have expected ensembl_gene_id/mgi_symbol columns; continuing with org.Mm.eg.db mapping.")
  }
} else {
  log_msg("No project gene_symbol_mapping.tsv found; continuing with org.Mm.eg.db mapping only.")
}

map_genes_to_entrez <- function(gene_ids, original_ids = gene_ids) {
  base_ids <- strip_ensembl_version(gene_ids)
  if (!has_pkg("org.Mm.eg.db")) {
    return(tibble(
      original_gene_id = original_ids,
      ensembl_base = base_ids,
      org_symbol = NA_character_,
      entrez_ids = NA_character_,
      chosen_entrez = NA_character_,
      mapped_flag = FALSE
    ))
  }
  raw_map <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = unique(base_ids),
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ENSEMBL"
  ) %>%
    as_tibble() %>%
    filter(!is.na(.data$ENSEMBL))
  summarised <- raw_map %>%
    group_by(ENSEMBL) %>%
    summarise(
      org_symbol = paste(sort(unique(na.omit(.data$SYMBOL))), collapse = ";"),
      entrez_ids = paste(sort(unique(na.omit(.data$ENTREZID))), collapse = ";"),
      chosen_entrez = dplyr::first(sort(unique(na.omit(.data$ENTREZID)))),
      .groups = "drop"
    ) %>%
    dplyr::rename(ensembl_base = ENSEMBL)
  tibble(original_gene_id = original_ids, ensembl_base = base_ids) %>%
    left_join(summarised, by = "ensembl_base") %>%
    mutate(
      org_symbol = na_if(.data$org_symbol, ""),
      entrez_ids = na_if(.data$entrez_ids, ""),
      chosen_entrez = as.character(.data$chosen_entrez),
      mapped_flag = !is.na(.data$chosen_entrez) & nzchar(.data$chosen_entrez)
    )
}

query_mapping <- map_genes_to_entrez(query_tbl$ensembl_base, query_tbl$original_gene_id) %>%
  left_join(project_map, by = "ensembl_base") %>%
  left_join(
    query_tbl %>%
      dplyr::select(query_order, original_gene_id, dplyr::any_of(c("weight", "beta_meta", "se_meta", "I2", "power"))),
    by = "original_gene_id"
  ) %>%
  relocate(query_order, original_gene_id, ensembl_base, project_mgi_symbol, org_symbol, entrez_ids, chosen_entrez, mapped_flag)

query_entrez <- query_mapping %>%
  filter(.data$mapped_flag) %>%
  pull(.data$chosen_entrez) %>%
  unique() %>%
  as.character()
n_query_mapped <- length(query_entrez)
n_query_unmapped <- sum(!query_mapping$mapped_flag)
log_msg("Mapped query Entrez IDs: ", n_query_mapped)
log_msg("Unmapped retained genes: ", n_query_unmapped)

write_tsv(query_mapping, file.path(tables_dir, "IMRS_retained_gene_mapping_audit.tsv"))

background_candidates <- list(
  gene_power = list(
    path = imrs_config_field_path(config, "gene_power"),
    definition = "All anchor-phase genes present in gene_power.tsv after anchor meta-analysis/power evaluation; used as the project-specific anchor-eligible background."
  ),
  gene_heterogeneity = list(
    path = imrs_config_field_path(config, "gene_heterogeneity"),
    definition = "All anchor-phase genes present in gene_heterogeneity.tsv after anchor heterogeneity evaluation; fallback project-specific background."
  ),
  support_by_dataset = list(
    path = imrs_config_field_path(config, "support_by_dataset"),
    definition = "All genes present in support_by_dataset.tsv across locked-anchor datasets; fallback anchor support background."
  ),
  core_support_all_genes = list(
    path = imrs_config_field_path(config, "core_support_all_genes"),
    definition = "All genes present in core_support_all_genes.tsv; fallback anchor support background."
  )
)

background_source_file <- NA_character_
background_definition <- NA_character_
background_original <- character()
for (nm in names(background_candidates)) {
  candidate <- background_candidates[[nm]]
  if (!file.exists(candidate$path)) next
  candidate_tbl <- readr::read_tsv(candidate$path, show_col_types = FALSE, progress = FALSE)
  bg_col <- detect_gene_column(candidate_tbl, c("gene_id", "gene", "ensembl_gene_id", "gene_id_clean", "ENSEMBL"))
  if ("phase" %in% names(candidate_tbl)) {
    candidate_tbl <- candidate_tbl %>% filter(is.na(.data$phase) | .data$phase == "anchor")
  }
  background_original <- candidate_tbl %>%
    transmute(gene_id = as.character(.data[[bg_col]])) %>%
    filter(!is.na(.data$gene_id), nzchar(.data$gene_id)) %>%
    distinct(.data$gene_id) %>%
    pull(.data$gene_id)
  if (length(background_original) > 0) {
    background_source_file <- candidate$path
    background_definition <- candidate$definition
    break
  }
}

if (length(background_original) == 0) {
  if (!has_pkg("org.Mm.eg.db")) {
    fail_with_log("No project-specific background found and org.Mm.eg.db unavailable for annotation fallback.")
  }
  log_msg("No project-specific background found; using weak org.Mm.eg.db annotation fallback.")
  background_source_file <- "org.Mm.eg.db"
  background_definition <- "Weak fallback: all mouse Entrez IDs available in org.Mm.eg.db."
  background_mapping <- tibble(
    original_gene_id = AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "ENTREZID"),
    ensembl_base = NA_character_,
    org_symbol = NA_character_,
    entrez_ids = original_gene_id,
    chosen_entrez = original_gene_id,
    mapped_flag = TRUE
  )
} else {
  log_msg("Background source: ", background_source_file)
  background_mapping <- map_genes_to_entrez(background_original, background_original)
}

background_entrez <- background_mapping %>%
  filter(.data$mapped_flag) %>%
  pull(.data$chosen_entrez) %>%
  unique() %>%
  as.character()
query_entrez_in_background <- intersect(query_entrez, background_entrez)
if (length(query_entrez_in_background) < length(query_entrez)) {
  log_msg("Mapped query IDs outside background universe: ", length(setdiff(query_entrez, background_entrez)))
}

background_audit <- tibble(
  background_source_file = background_source_file,
  background_definition = background_definition,
  n_background_original = length(unique(background_original)),
  n_background_mapped_entrez = length(background_entrez),
  n_query_original = n_query_original,
  n_query_mapped_entrez = n_query_mapped,
  n_query_unmapped = n_query_unmapped,
  notes = paste0(
    "Query genes in mapped background: ", length(query_entrez_in_background),
    ". Enrichment universe uses mapped background Entrez IDs. No IMRS scoring or retained-gene definitions were changed."
  )
)
write_tsv(background_audit, file.path(tables_dir, "IMRS_enrichment_background_audit.tsv"))

analysis_background_label <- if (grepl("gene_power.tsv$", background_source_file)) {
  "project-specific anchor-eligible gene_power.tsv background"
} else if (grepl("gene_heterogeneity.tsv$", background_source_file)) {
  "project-specific anchor gene_heterogeneity.tsv background"
} else if (grepl("support_by_dataset.tsv$", background_source_file)) {
  "project-specific anchor support_by_dataset.tsv background"
} else {
  background_definition
}

symbol_by_entrez <- query_mapping %>%
  filter(.data$mapped_flag) %>%
  transmute(chosen_entrez = as.character(.data$chosen_entrez),
            symbol = dplyr::coalesce(na_if(.data$org_symbol, ""), na_if(.data$project_mgi_symbol, ""), .data$ensembl_base)) %>%
  distinct(.data$chosen_entrez, .keep_all = TRUE)

symbols_for_entrez <- function(entrez_string) {
  ids <- unlist(strsplit(as.character(entrez_string), "/", fixed = TRUE))
  ids <- ids[nzchar(ids)]
  syms <- symbol_by_entrez$symbol[match(ids, symbol_by_entrez$chosen_entrez)]
  syms <- ifelse(is.na(syms) | !nzchar(syms), ids, syms)
  paste(unique(syms), collapse = ";")
}

parse_ratio_num <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  as.integer(vapply(parts, function(z) if (length(z) >= 1) z[[1]] else NA_character_, character(1)))
}

format_enrichment_table <- function(raw_df, readable_df, database, package_name) {
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    return(tibble(
      database = character(),
      term_id = character(),
      term_name = character(),
      description = character(),
      gene_count = integer(),
      background_count = integer(),
      gene_ratio = character(),
      bg_ratio = character(),
      p_value = numeric(),
      adjusted_p_value = numeric(),
      q_value = numeric(),
      overlapping_genes_symbols = character(),
      overlapping_genes_entrez = character(),
      analysis_background = character(),
      query_gene_count = integer(),
      mapped_query_gene_count = integer(),
      background_gene_count = integer(),
      software_package = character(),
      software_version = character()
    ))
  }
  raw_tbl <- as_tibble(raw_df)
  readable_tbl <- if (!is.null(readable_df) && nrow(readable_df) > 0 && "geneID" %in% names(readable_df)) {
    as_tibble(readable_df) %>% dplyr::select(ID, readable_geneID = geneID)
  } else {
    tibble(ID = raw_tbl$ID, readable_geneID = vapply(raw_tbl$geneID, symbols_for_entrez, character(1)))
  }
  raw_tbl %>%
    left_join(readable_tbl, by = "ID") %>%
    transmute(
      database = database,
      term_id = .data$ID,
      term_name = .data$Description,
      description = .data$Description,
      gene_count = as.integer(.data$Count),
      background_count = parse_ratio_num(.data$BgRatio),
      gene_ratio = .data$GeneRatio,
      bg_ratio = .data$BgRatio,
      p_value = .data$pvalue,
      adjusted_p_value = .data$p.adjust,
      q_value = if ("qvalue" %in% names(.)) .data$qvalue else NA_real_,
      overlapping_genes_symbols = stringr::str_replace_all(dplyr::coalesce(.data$readable_geneID, ""), "/", ";"),
      overlapping_genes_entrez = stringr::str_replace_all(.data$geneID, "/", ";"),
      analysis_background = analysis_background_label,
      query_gene_count = n_query_original,
      mapped_query_gene_count = n_query_mapped,
      background_gene_count = length(background_entrez),
      software_package = package_name,
      software_version = as.character(utils::packageVersion(package_name))
    ) %>%
    arrange(.data$adjusted_p_value, .data$p_value, dplyr::desc(.data$gene_count))
}

run_go <- function() {
  if (!all(vapply(c("clusterProfiler", "org.Mm.eg.db"), has_pkg, logical(1)))) {
    log_msg("GO BP unavailable; missing package(s): ",
            paste(c("clusterProfiler", "org.Mm.eg.db")[!vapply(c("clusterProfiler", "org.Mm.eg.db"), has_pkg, logical(1))], collapse = ", "))
    return(list(table = format_enrichment_table(NULL, NULL, "GO Biological Process", "clusterProfiler"), object = NULL))
  }
  log_msg("Running GO Biological Process over-representation analysis.")
  raw <- tryCatch(
    clusterProfiler::enrichGO(
      gene = query_entrez_in_background,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      ont = "BP",
      keyType = "ENTREZID",
      universe = background_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    ),
    error = function(e) {
      log_msg("GO BP analysis failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw)) return(list(table = format_enrichment_table(NULL, NULL, "GO Biological Process", "clusterProfiler"), object = NULL))
  readable <- tryCatch(clusterProfiler::setReadable(raw, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID"),
                       error = function(e) NULL)
  raw_df <- as.data.frame(raw)
  readable_df <- if (is.null(readable)) NULL else as.data.frame(readable)
  list(table = format_enrichment_table(raw_df, readable_df, "GO Biological Process", "clusterProfiler"), object = raw)
}

run_reactome <- function() {
  if (!all(vapply(c("ReactomePA"), has_pkg, logical(1)))) {
    log_msg("Reactome unavailable; missing package: ReactomePA")
    return(list(table = format_enrichment_table(NULL, NULL, "Reactome", "ReactomePA"), object = NULL))
  }
  log_msg("Running Reactome pathway over-representation analysis.")
  raw <- tryCatch(
    ReactomePA::enrichPathway(
      gene = query_entrez_in_background,
      organism = "mouse",
      universe = background_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    ),
    error = function(e) {
      log_msg("Reactome analysis failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw)) return(list(table = format_enrichment_table(NULL, NULL, "Reactome", "ReactomePA"), object = NULL))
  readable <- tryCatch(clusterProfiler::setReadable(raw, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID"),
                       error = function(e) NULL)
  raw_df <- as.data.frame(raw)
  readable_df <- if (is.null(readable)) NULL else as.data.frame(readable)
  list(table = format_enrichment_table(raw_df, readable_df, "Reactome", "ReactomePA"), object = raw)
}

load_hallmark_term2gene <- function() {
  if (!has_pkg("msigdbr")) return(NULL)
  collections <- tryCatch(msigdbr::msigdbr_collections(db_species = "MM"), error = function(e) NULL)
  collection_arg <- if (!is.null(collections) && "MH" %in% collections$gs_collection) "MH" else "H"
  hallmark <- tryCatch(
    msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = collection_arg),
    error = function(e) {
      log_msg("Mouse-native MSigDB Hallmark load failed with collection ", collection_arg, ": ", conditionMessage(e))
      tryCatch(msigdbr::msigdbr(db_species = "HS", species = "Mus musculus", collection = "H"),
               error = function(e2) {
                 log_msg("Human Hallmark ortholog fallback load failed: ", conditionMessage(e2))
                 NULL
               })
    }
  )
  if (is.null(hallmark) || nrow(hallmark) == 0) return(NULL)
  id_col <- if ("ncbi_gene" %in% names(hallmark)) "ncbi_gene" else if ("entrez_gene" %in% names(hallmark)) "entrez_gene" else NA_character_
  if (is.na(id_col)) {
    log_msg("MSigDB Hallmark table did not include an Entrez/NCBI gene column.")
    return(NULL)
  }
  hallmark %>%
    transmute(term = .data$gs_name, gene = as.character(.data[[id_col]])) %>%
    filter(!is.na(.data$term), !is.na(.data$gene), nzchar(.data$gene)) %>%
    distinct()
}

run_hallmark <- function() {
  if (!all(vapply(c("clusterProfiler", "msigdbr"), has_pkg, logical(1)))) {
    log_msg("MSigDB Hallmark unavailable; missing package(s): ",
            paste(c("clusterProfiler", "msigdbr")[!vapply(c("clusterProfiler", "msigdbr"), has_pkg, logical(1))], collapse = ", "))
    return(list(table = format_enrichment_table(NULL, NULL, "MSigDB Hallmark", "clusterProfiler"), object = NULL))
  }
  term2gene <- load_hallmark_term2gene()
  if (is.null(term2gene) || nrow(term2gene) == 0) {
    log_msg("MSigDB Hallmark unavailable; no TERM2GENE table could be loaded.")
    return(list(table = format_enrichment_table(NULL, NULL, "MSigDB Hallmark", "clusterProfiler"), object = NULL))
  }
  log_msg("Running MSigDB Hallmark over-representation analysis with ", n_distinct(term2gene$term), " Hallmark gene sets.")
  hallmark_universe <- intersect(background_entrez, unique(term2gene$gene))
  hallmark_query <- intersect(query_entrez_in_background, hallmark_universe)
  raw <- tryCatch(
    clusterProfiler::enricher(
      gene = hallmark_query,
      universe = hallmark_universe,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 10,
      maxGSSize = 500
    ),
    error = function(e) {
      log_msg("MSigDB Hallmark analysis failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw)) return(list(table = format_enrichment_table(NULL, NULL, "MSigDB Hallmark", "clusterProfiler"), object = NULL))
  raw_df <- as.data.frame(raw)
  list(table = format_enrichment_table(raw_df, NULL, "MSigDB Hallmark", "clusterProfiler"), object = raw)
}

go_result <- run_go()
reactome_result <- run_reactome()
hallmark_result <- run_hallmark()

all_enrichment <- bind_rows(go_result$table, reactome_result$table, hallmark_result$table)
s5_tsv <- file.path(tables_dir, "Supplementary_Table_S5_IMRS_gene_enrichment_all.tsv")
write_tsv(all_enrichment, s5_tsv)
if (has_pkg("writexl")) {
  writexl::write_xlsx(
    list(
      enrichment_all = all_enrichment,
      retained_gene_mapping_audit = query_mapping,
      background_audit = background_audit
    ),
    file.path(tables_dir, "Supplementary_Table_S5_IMRS_gene_enrichment_all.xlsx")
  )
  log_msg("Wrote Supplementary Table S5 TSV and XLSX.")
} else {
  log_msg("writexl/openxlsx unavailable; wrote TSV only. Install command: ", install_commands[["writexl"]])
}

sig_count <- function(db) {
  all_enrichment %>% filter(.data$database == db, !is.na(.data$adjusted_p_value), .data$adjusted_p_value <= 0.05) %>% nrow()
}
go_n <- sig_count("GO Biological Process")
reactome_n <- sig_count("Reactome")
hallmark_n <- sig_count("MSigDB Hallmark")
log_msg("Significant GO BP terms (FDR <= 0.05): ", go_n)
log_msg("Significant Reactome terms (FDR <= 0.05): ", reactome_n)
log_msg("Significant MSigDB Hallmark terms (FDR <= 0.05): ", hallmark_n)

theme_imrs_enrichment <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#111827", margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(face = "bold", color = "#111827"),
      axis.text = ggplot2::element_text(color = "#1F2937"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 16, 10, 10)
    )
}

plot_enrichment <- function(tbl, database_label, title, placeholder_label) {
  hallmark_readable_label <- function(term_name) {
    cleaned <- stringr::str_remove(term_name, "^HALLMARK_")
    dplyr::case_when(
      term_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE" ~ "Interferon gamma response",
      term_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE" ~ "Interferon alpha response",
      term_name == "HALLMARK_INFLAMMATORY_RESPONSE" ~ "Inflammatory response",
      term_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING" ~ "IL6/JAK/STAT3 signaling",
      term_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~ "TNFα/NF-κB signaling",
      term_name == "HALLMARK_ALLOGRAFT_REJECTION" ~ "Allograft rejection",
      term_name == "HALLMARK_APOPTOSIS" ~ "Apoptosis",
      TRUE ~ stringr::str_to_sentence(stringr::str_replace_all(cleaned, "_", " "))
    )
  }
  plot_tbl <- tbl %>%
    filter(.data$database == .env$database_label, !is.na(.data$adjusted_p_value), .data$adjusted_p_value <= 0.05) %>%
    arrange(.data$adjusted_p_value, .data$p_value, desc(.data$gene_count)) %>%
    slice_head(n = 15) %>%
    mutate(
      plot_term_name = if (identical(.env$database_label, "MSigDB Hallmark")) hallmark_readable_label(.data$term_name) else .data$term_name,
      term_wrapped = stringr::str_wrap(.data$plot_term_name, width = 42),
      neg_log10_fdr = -log10(pmax(.data$adjusted_p_value, .Machine$double.xmin)),
      term_wrapped = forcats::fct_reorder(.data$term_wrapped, .data$neg_log10_fdr)
    )
  if (nrow(plot_tbl) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = placeholder_label, size = 4.0, color = "#374151") +
        labs(title = title) +
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
        theme_void(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", color = "#111827", margin = margin(b = 8)),
          plot.margin = margin(10, 16, 10, 10)
        )
    )
  }
  ggplot(plot_tbl, aes(x = .data$gene_count, y = .data$term_wrapped)) +
    geom_point(aes(size = .data$gene_count, color = .data$neg_log10_fdr), alpha = 0.92) +
    scale_color_gradient(low = "#7AA6C2", high = "#123B5D", name = "-log10(FDR)") +
    scale_size_continuous(name = "Gene count", range = c(2.5, 7)) +
    labs(title = title, x = "Gene count", y = NULL) +
    theme_imrs_enrichment()
}

save_plot_all_formats <- function(plot, stem, width, height) {
  png_path <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figures_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(figures_dir, paste0(stem, ".svg"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 400, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  if (has_pkg("svglite")) {
    ggplot2::ggsave(svg_path, plot, width = width, height = height, device = svglite::svglite, bg = "white", limitsize = FALSE)
  } else {
    log_msg("svglite unavailable; SVG skipped for ", stem, ". Install command: ", install_commands[["svglite"]])
  }
  c(png = png_path, pdf = pdf_path, svg = svg_path)
}

plot_go <- plot_enrichment(
  all_enrichment,
  "GO Biological Process",
  "GO Biological Process enrichment of retained frozen IMRS genes",
  "No terms passed FDR threshold for GO Biological Process."
)
plot_reactome <- plot_enrichment(
  all_enrichment,
  "Reactome",
  "Reactome enrichment of retained frozen IMRS genes",
  "No terms passed FDR threshold for Reactome."
)
plot_hallmark <- plot_enrichment(
  all_enrichment,
  "MSigDB Hallmark",
  "MSigDB Hallmark enrichment of retained frozen IMRS genes",
  "No terms passed FDR threshold for MSigDB Hallmark."
)

paths_go <- save_plot_all_formats(plot_go, "FigureS2A_GO_BP_enrichment", 7.5, 5.4)
paths_reactome <- save_plot_all_formats(plot_reactome, "FigureS2B_Reactome_enrichment", 7.5, 5.4)
paths_hallmark <- save_plot_all_formats(plot_hallmark, "FigureS2C_MSigDB_Hallmark_enrichment", 7.5, 5.4)
combined_plot <- (plot_go / plot_reactome / plot_hallmark) +
  patchwork::plot_layout(heights = c(1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A")
paths_combined <- save_plot_all_formats(combined_plot, "FigureS2_gene_program_enrichment_combined", 8.8, 14.8)
log_msg("Wrote Figure S2 panel and combined figure outputs.")

caption_text <- paste0(
  "Supplementary Figure S2. Gene-program enrichment of retained frozen IMRS genes supports an acute delivery-associated innate transcriptional program. ",
  "(A) GO Biological Process over-representation analysis of the 287 retained frozen IMRS genes. ",
  "(B) Reactome pathway enrichment analysis of the same retained gene set. ",
  "(C) MSigDB Hallmark enrichment analysis. ",
  "Dot size indicates the number of retained IMRS genes overlapping each term, and color indicates −log10(FDR). ",
  "Enriched terms were concentrated in antiviral, interferon-associated, cytokine/chemokine, inflammatory, and innate immune response programs. ",
  "These analyses provide program-level biological context for the frozen IMRS coefficients but do not establish causal pathways, cell-type sources, clinical reactogenicity prediction, or delivery-platform safety ranking."
)
writeLines(caption_text, file.path(analysis_root, "FigureS2_caption_candidate.txt"), useBytes = TRUE)

results_paragraph <- paste0(
  "To further evaluate the biological structure of the retained IMRS gene set, we performed gene-program enrichment analysis using the 287 frozen IMRS genes against an anchor-eligible background gene universe. ",
  "Enriched terms were concentrated in antiviral response, interferon-associated signaling, cytokine/chemokine signaling, inflammatory response, and innate immune response programs (Supplementary Figure S2; Supplementary Table S5). ",
  "These results support interpreting the frozen coefficients as an acute delivery-associated innate transcriptional program, while remaining consistent with the scope of IMRS as a transfer-oriented scoring framework rather than a mechanistic pathway model."
)
mapping_note <- "Mapping note for Methods/Supplement: Of 287 retained IMRS genes, 252 mapped to Entrez identifiers and 35 were unmapped. The enrichment background was derived from the anchor gene_power.tsv universe, with 22,304 mapped background genes."
writeLines(c(results_paragraph, "", mapping_note), file.path(analysis_root, "Results_insert_candidate.txt"), useBytes = TRUE)

file_size_summary <- function(paths) {
  labels <- names(paths)
  paste(vapply(seq_along(paths), function(i) {
    p <- paths[[i]]
    paste0(labels[[i]], "=", if (file.exists(p)) file.info(p)$size else 0)
  }, character(1)), collapse = ";")
}

file_exists_summary <- function(paths) {
  labels <- names(paths)
  paste(vapply(seq_along(paths), function(i) {
    p <- paths[[i]]
    paste0(labels[[i]], "=", file.exists(p) && file.info(p)$size > 0)
  }, character(1)), collapse = ";")
}

append_or_replace_change_log <- function() {
  path <- file.path(v6_root, "v6_change_log.tsv")
  cols <- c("file", "panel_or_figure", "change_type", "old_text", "new_text",
            "plot_visible_or_caption", "regenerated_or_copied", "notes")
  existing <- if (file.exists(path)) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) else tibble::tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>%
    filter(!(.data$file == "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R" &
               .data$change_type == "analysis_added")) %>%
    filter(.data$change_type != "Priority3_enrichment_presentation_update")
  analysis_row <- tibble(
    file = "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R",
    panel_or_figure = "Supplementary Figure S2 / Supplementary Table S5",
    change_type = "analysis_added",
    old_text = "",
    new_text = "Added Priority 3 retained-gene enrichment analysis of 287 frozen IMRS genes.",
    plot_visible_or_caption = "analysis/figure/table",
    regenerated_or_copied = "generated",
    notes = "No IMRS scoring framework, retained genes, validation, robustness, or transfer-evaluation numbers were changed."
  )
  presentation_rows <- tibble::tribble(
    ~file, ~panel_or_figure, ~change_type, ~old_text, ~new_text, ~plot_visible_or_caption, ~regenerated_or_copied, ~notes,
    "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R", "Supplementary Figure S2A", "Priority3_enrichment_presentation_update", "GO Biological Process enrichment of retained IMRS genes", "GO Biological Process enrichment of retained frozen IMRS genes", "plot_visible", "regenerated", "Presentation-only title update; enrichment statistics unchanged.",
    "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R", "Supplementary Figure S2B", "Priority3_enrichment_presentation_update", "Reactome enrichment of retained IMRS genes", "Reactome enrichment of retained frozen IMRS genes", "plot_visible", "regenerated", "Presentation-only title update; enrichment statistics unchanged.",
    "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R", "Supplementary Figure S2C", "Priority3_enrichment_presentation_update", "MSigDB Hallmark enrichment of retained IMRS genes", "MSigDB Hallmark enrichment of retained frozen IMRS genes", "plot_visible", "regenerated", "Presentation-only title update; enrichment statistics unchanged.",
    "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R", "Supplementary Figure S2C", "Priority3_enrichment_presentation_update", "Raw all-caps MSigDB Hallmark IDs in plot labels", "Readable Hallmark plot labels; exact MSigDB term IDs retained in Supplementary Table S5", "plot_visible", "regenerated", "Figure label readability update only; enrichment statistics unchanged.",
    "Priority3_gene_program_enrichment/FigureS2_caption_candidate.txt", "Supplementary Figure S2", "Priority3_enrichment_presentation_update", "Previous Supplementary Figure S2 caption candidate", "Revised caption naming retained frozen IMRS genes and explaining dot size/color.", "caption", "generated", "Caption candidate updated for manuscript-review suggestions.",
    "Priority3_gene_program_enrichment/Results_insert_candidate.txt", "Supplementary Figure S2 / Supplementary Table S5", "Priority3_enrichment_presentation_update", "Previous Results insert candidate", "Revised Results paragraph plus separate mapping note.", "caption/text", "generated", "Mapping note kept below the Results paragraph.",
    "Priority3_gene_program_enrichment/Priority3_enrichment_analysis_log.txt", "Supplementary Figure S2", "Priority3_enrichment_presentation_update", "", "Logged presentation-only title and Hallmark readability updates.", "log", "generated", "No enrichment statistics, input genes, background universe, IMRS scoring, or validation results changed."
  )
  write_tsv(bind_rows(existing[, cols], analysis_row, presentation_rows), path)
}

append_or_replace_output_manifest <- function() {
  path <- file.path(v6_root, "v6_output_manifest.tsv")
  cols <- c("figure_id", "panel_id", "output_png", "output_pdf", "output_svg",
            "source_script", "regenerated_or_copied", "file_exists", "file_size_bytes", "notes")
  existing <- if (file.exists(path)) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) else tibble::tibble()
  for (col in cols) if (!col %in% names(existing)) existing[[col]] <- character()
  existing <- existing %>%
    filter(!(.data$figure_id == "Supplementary Figure S2"))
  mkrow <- function(panel_id, paths, notes) {
    tibble(
      figure_id = "Supplementary Figure S2",
      panel_id = panel_id,
      output_png = normalizePath(paths[["png"]], winslash = "/", mustWork = FALSE),
      output_pdf = normalizePath(paths[["pdf"]], winslash = "/", mustWork = FALSE),
      output_svg = normalizePath(paths[["svg"]], winslash = "/", mustWork = FALSE),
      source_script = "Priority3_gene_program_enrichment/run_IMRS_retained_gene_enrichment_v6.R",
      regenerated_or_copied = "generated",
      file_exists = file_exists_summary(paths),
      file_size_bytes = file_size_summary(paths),
      notes = notes
    )
  }
  rows <- bind_rows(
    mkrow("A", paths_go, "GO Biological Process enrichment of retained frozen IMRS genes."),
    mkrow("B", paths_reactome, "Reactome enrichment of retained frozen IMRS genes."),
    mkrow("C", paths_hallmark, "MSigDB Hallmark enrichment of retained frozen IMRS genes with readable plot labels; exact MSigDB IDs retained in Supplementary Table S5."),
    mkrow("combined", paths_combined, "Combined Supplementary Figure S2 gene-program enrichment panels.")
  )
  write_tsv(bind_rows(existing[, cols], rows), path)
}

infer_table_purpose <- function(rel) {
  dplyr::case_when(
    grepl("Supplementary_Table_S5_IMRS_gene_enrichment_all\\.tsv$", rel) ~ "Priority 3 enrichment result table",
    grepl("Supplementary_Table_S5_IMRS_gene_enrichment_all\\.xlsx$", rel) ~ "Priority 3 enrichment result workbook",
    grepl("IMRS_retained_gene_mapping_audit\\.tsv$", rel) ~ "Priority 3 retained-gene ID mapping audit",
    grepl("IMRS_enrichment_background_audit\\.tsv$", rel) ~ "Priority 3 enrichment background audit",
    TRUE ~ "v6 figure-supporting table"
  )
}

update_table_inventory <- function() {
  path <- file.path(v6_root, "v6_table_inventory.tsv")
  files <- c(
    list.files(v6_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE),
    list.files(v6_root, pattern = "Supplementary_Table_S5_IMRS_gene_enrichment_all\\.xlsx$", recursive = TRUE, full.names = TRUE)
  )
  rows <- purrr::map_dfr(sort(files), function(p) {
    ext <- tolower(tools::file_ext(p))
    if (identical(ext, "tsv")) {
      lines <- readLines(p, warn = FALSE, encoding = "UTF-8")
      row_count <- max(length(lines) - 1L, 0L)
      column_count <- if (length(lines) > 0) length(strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]) else 0L
    } else if (grepl("Supplementary_Table_S5_IMRS_gene_enrichment_all\\.xlsx$", p)) {
      row_count <- nrow(all_enrichment)
      column_count <- ncol(all_enrichment)
    } else {
      row_count <- NA_integer_
      column_count <- NA_integer_
    }
    rel <- stringr::str_replace(normalizePath(p, winslash = "/", mustWork = FALSE),
                                paste0("^", stringr::fixed(normalizePath(v6_root, winslash = "/", mustWork = FALSE)), "/?"), "")
    source_path <- file.path(imrs_config_field_path(config, "source_v5_dir"), rel)
    tibble(
      table_file = rel,
      purpose_inferred = infer_table_purpose(rel),
      edited_yes_no = ifelse(!file.exists(source_path), "yes", "no"),
      key_changes = ifelse(grepl("^Priority3_gene_program_enrichment/", rel), "Created for Priority 3 retained-gene enrichment analysis", ""),
      row_count = row_count,
      column_count = column_count
    )
  })
  write_tsv(rows, path)
}

update_file_inventory <- function() {
  path <- file.path(v6_root, "v6_file_inventory.tsv")
  files <- list.files(v6_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  rows <- purrr::map_dfr(sort(files), function(p) {
    info <- file.info(p)
    rel <- stringr::str_replace(normalizePath(p, winslash = "/", mustWork = FALSE),
                                paste0("^", stringr::fixed(normalizePath(v6_root, winslash = "/", mustWork = FALSE)), "/?"), "")
    source_path <- file.path(imrs_config_field_path(config, "source_v5_dir"), rel)
    same_as_v5 <- file.exists(source_path) && identical(unname(tools::md5sum(p)), unname(tools::md5sum(source_path)))
    tibble(
      file = rel,
      extension = tools::file_ext(p),
      size_bytes = as.numeric(info$size),
      modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
      copied_from_v5_yes_no = ifelse(file.exists(source_path), "yes", "no"),
      generated_in_v6_yes_no = ifelse(!same_as_v5, "yes", "no")
    )
  })
  write_tsv(rows, path)
}

append_or_replace_change_log()
append_or_replace_output_manifest()
update_table_inventory()
update_file_inventory()
log_msg("Updated v6 change log, output manifest, table inventory, and file inventory.")
log_msg("Presentation update: plot titles were updated to 'retained frozen IMRS genes'.")
log_msg("Presentation update: Hallmark figure labels were converted to readable labels while exact MSigDB IDs remain preserved in Supplementary Table S5.")
log_msg("Presentation update: no enrichment statistics, input genes, background universe, IMRS scoring, or validation results were changed.")

validation_lines <- c(
  paste0("Query retained gene count exactly 287: ", n_query_original == 287L),
  paste0("Query mapped Entrez gene count reported: ", n_query_mapped),
  paste0("Background source and mapped size reported: ", background_source_file, " / ", length(background_entrez)),
  paste0("GO BP generated or logged unavailable: TRUE; significant terms=", go_n),
  paste0("Reactome generated or logged unavailable: TRUE; significant terms=", reactome_n),
  paste0("MSigDB Hallmark generated or logged unavailable: TRUE; significant terms=", hallmark_n),
  paste0("Supplementary Table S5 TSV exists: ", file.exists(s5_tsv)),
  paste0("Supplementary Table S5 XLSX exists: ", file.exists(file.path(tables_dir, "Supplementary_Table_S5_IMRS_gene_enrichment_all.xlsx"))),
  paste0("All generated figure files nonzero: ", all(file.exists(c(paths_go, paths_reactome, paths_hallmark, paths_combined)) &
                                                       file.info(c(paths_go, paths_reactome, paths_hallmark, paths_combined))$size > 0)),
  "No IMRS score, validation result, robustness number, transfer-evaluation number, or existing figure data was modified by this script.",
  "No claims of mechanism, cell-type source, clinical prediction, or delivery-platform safety ranking were introduced."
)
log_msg("Validation summary:")
for (line in validation_lines) log_msg("  ", line)
log_msg("Priority 3 retained IMRS gene-program enrichment analysis complete.")
write_log()

cat("Priority 3 enrichment complete\n")
cat("Retained query genes: ", n_query_original, "\n", sep = "")
cat("Mapped query Entrez genes: ", n_query_mapped, "\n", sep = "")
cat("Background: ", analysis_background_label, " (mapped Entrez n=", length(background_entrez), ")\n", sep = "")
cat("Significant GO BP terms: ", go_n, "\n", sep = "")
cat("Significant Reactome terms: ", reactome_n, "\n", sep = "")
cat("Significant MSigDB Hallmark terms: ", hallmark_n, "\n", sep = "")
