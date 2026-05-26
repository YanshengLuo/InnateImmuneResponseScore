#!/usr/bin/env Rscript

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "08_manuscript_ready_cleanup"

claim_group_map <- c(
  anchor = "anchor_discovery",
  primary_acute_validation = "primary_external_validation",
  extended_validation = "extended_exploratory_transfer",
  secondary_support = "secondary_support_not_primary",
  excluded_or_unclear = "unclear_not_for_main_claim"
)

risk_levels <- c(high = 3L, medium = 2L, low = 1L)

locked_anchor_dataset_ids <- c(
  "GSE39129", "GSE167521", "GSE264344", "GSE279372", "GSE279744"
)
production_discovery_datasets <- locked_anchor_dataset_ids

role_display_map <- c(
  anchor = "Locked anchor",
  primary_acute_validation = "Primary acute validation",
  extended_validation = "Extended validation",
  secondary_support = "Secondary support",
  excluded_or_unclear = "Excluded/unclear"
)

base_gse_id <- function(dataset_id) {
  out <- str_match(as.character(dataset_id), "^(GSE[0-9]+)")[, 2]
  ifelse(is.na(out), as.character(dataset_id), out)
}

manual_metadata <- function() {
  tibble::tribble(
    ~dataset_id, ~gse_id_clean, ~organism_clean, ~delivery_platform_clean, ~PMID_clean, ~paper_title_clean, ~DOI_clean, ~general_biology, ~role_note_manual,
    "GSE39129", "GSE39129", "Mus musculus", "lentiviral vector", NA_character_, NA_character_, NA_character_, "acute systemic viral vector delivery", "locked anchor; mouse liver; 4 h",
    "GSE167521", "GSE167521", "Mus musculus", "lipid nanoparticle", NA_character_, NA_character_, NA_character_, "acute local LNP delivery", "locked anchor; mouse skin; 24 h",
    "GSE264344", "GSE264344", "Mus musculus", "adenoviral vector", "40162786", "Early spatiotemporal evolution of the immune response elicited by adenovirus serotype 26 vector vaccination in mice", "10.1128/jvi.00247-25", "adenoviral vector vaccination time-course across blood, muscle, and draining lymph node", NA_character_,
    "GSE279744", "GSE279744", "Mus musculus", "mRNA lipid nanoparticle vaccine", NA_character_, NA_character_, NA_character_, "LNP formulation calibration; lymph node; 6 h", "locked anchor; acute LNP formulation / dLN context",
    "GSE279372", "GSE279372", "Mus musculus", "vaccine formulation / LNP or alum/VLP depending split", NA_character_, NA_character_, NA_character_, "acute draining lymph node formulation response", "locked anchor; acute discovery context",
    "GSE262515_cell_line", "GSE262515", "Homo sapiens", "optimized therapeutic lipid nanoparticle", "38902241", "Optimized RNA interference therapeutics combined with interleukin-2 mRNA for treating hepatitis B virus infection", "10.1038/s41392-024-01871-8", "optimized tLNP therapeutic RNAi context; cell-line arm; 72 h", "secondary support / weak-context audit; not in production discovery set",
    "GSE262515_tissue", "GSE262515", "Mus musculus", "optimized therapeutic lipid nanoparticle", "38902241", "Optimized RNA interference therapeutics combined with interleukin-2 mRNA for treating hepatitis B virus infection", "10.1038/s41392-024-01871-8", "optimized tLNP therapeutic RNAi context; mouse tissue/liver arm; 72 h", "secondary support / weak-context audit; not in production discovery set",
    "GSE119119", "GSE119119", "Mus musculus", "adenoviral vector / virus-antibody perturbation", NA_character_, NA_character_, NA_character_, "acute liver viral perturbation; TRIM21-related context", "external acute validation, but interpret with genotype/antibody perturbation context",
    "GSE139529", "GSE139529", "Mus musculus", "adenoviral / MVA vaccine vector", NA_character_, NA_character_, NA_character_, "blood response 24 h after vaccination", "external acute validation",
    "GSE166655", "GSE166655", "Mus musculus", "AAV9 gene therapy", "34417184", "Gene therapy with AR isoform 2 rescues spinal and bulbar muscular atrophy phenotype by modulating AR transcriptional activity", "10.1126/sciadv.abi6896", "late AAV therapeutic transgene/disease-rescue model", "external extended transfer; not acute innate validation",
    "GSE178313", "GSE178313", "Mus musculus or humanized mouse if metadata indicates", "lipid nanoparticle", NA_character_, NA_character_, NA_character_, "liver LNP response at 48 h", "external extended transfer; not strict acute validation",
    "GSE279743", "GSE279743", "Mus musculus", "mRNA lipid nanoparticle vaccine", NA_character_, NA_character_, NA_character_, "lymph node response at 6 h", "external acute validation",
    "GSE314070", "GSE314070", "Mus musculus", "mRNA lipid nanoparticle vaccine", "41671339", "Parenteral vaccination with an adjuvanted mRNA vaccine induces protective mucosal immunity against rotavirus in neonatal mice", "10.1126/scitranslmed.adw6105", "late small-intestine mucosal/adaptive response after booster immunization", "external extended transfer; not acute innate validation"
  )
}

is_anchor_phase_row <- function(time_h, split_path = NA_character_) {
  time_h <- suppressWarnings(as.numeric(time_h))
  path_has_anchor <- str_detect(
    str_replace_all(as.character(split_path), "\\\\", "/"),
    "/deseq2_contrasts/anchor(/|$)"
  )
  (is.finite(time_h) & time_h <= 24) | replace_na(path_has_anchor, FALSE)
}

assign_manuscript_role <- function(dataset_id, time_h, split_path = NA_character_) {
  dataset_id <- as.character(dataset_id)
  time_h <- suppressWarnings(as.numeric(time_h))
  anchor_phase <- is_anchor_phase_row(time_h, split_path)
  case_when(
    dataset_id %in% locked_anchor_dataset_ids & anchor_phase ~ "anchor",
    dataset_id %in% c("GSE119119", "GSE139529", "GSE279743") ~ "primary_acute_validation",
    dataset_id == "GSE264344" & is.finite(time_h) & time_h > 24 ~ "extended_validation",
    dataset_id %in% c("GSE166655", "GSE178313", "GSE314070") ~ "extended_validation",
    dataset_id %in% c("GSE262515_cell_line", "GSE262515_tissue") ~ "secondary_support",
    is.finite(time_h) & time_h > 24 ~ "extended_validation",
    is.finite(time_h) & time_h <= 24 ~ "primary_acute_validation",
    TRUE ~ "excluded_or_unclear"
  )
}

role_note_for_row <- function(dataset_id, time_h, role_note_manual) {
  if (dataset_id == "GSE264344") {
    if (is.finite(time_h) && time_h <= 24) return("time_h <= 24: locked anchor discovery evidence")
    if (is.finite(time_h) && time_h > 24) {
      return("time_h > 24: extended validation / waning tissue-time kinetics")
    }
  }
  role_note_manual
}

safe_summary <- function(df) {
  tibble(
    n_contrasts = nrow(df),
    gse_ids = collapse_unique_chr(df$gse_id),
    mean_delta_imrs_z = safe_mean(df$delta_mean_imrs_z),
    median_delta_imrs_z = safe_median(df$delta_mean_imrs_z),
    proportion_positive_delta = safe_prop(df$delta_mean_imrs_z > 0),
    mean_auc_imrs_z_secondary = safe_mean(df$auc_imrs_z),
    median_auc_imrs_z_secondary = safe_median(df$auc_imrs_z)
  )
}

status_from_prop <- function(prop, pass_cutoff = 0.75, warn_cutoff = 0.50) {
  if (!is.finite(prop)) return("FAIL")
  if (prop >= pass_cutoff) return("PASS")
  if (prop >= warn_cutoff) return("WARNING")
  "FAIL"
}

risk_from_status <- function(status) {
  case_when(
    status == "PASS" ~ "low",
    status == "WARNING" ~ "medium",
    TRUE ~ "high"
  )
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

table_to_html <- function(tbl, title, file_name = NULL, max_rows = Inf) {
  shown <- tbl
  note <- ""
  if (is.finite(max_rows) && nrow(tbl) > max_rows) {
    shown <- head(tbl, max_rows)
    note <- paste0("<p><em>Showing first ", max_rows, " of ", nrow(tbl), " rows.</em></p>")
  }
  header <- paste0("<tr>", paste0("<th>", html_escape(names(shown)), "</th>", collapse = ""), "</tr>")
  rows <- if (nrow(shown) == 0) {
    "<tr><td>NA</td></tr>"
  } else {
    paste(apply(as.data.frame(shown), 1, function(r) {
      paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>")
    }), collapse = "\n")
  }
  link <- if (!is.null(file_name)) {
    paste0("<p><a href=\"../results/", html_escape(file_name), "\">",
           html_escape(file_name), "</a> (", nrow(tbl), " rows)</p>")
  } else {
    paste0("<p>", nrow(tbl), " rows</p>")
  }
  paste0("<h3>", html_escape(title), "</h3>", link, note,
         "<div class=\"table-wrap\"><table>", header, rows, "</table></div>")
}

plot_html <- function(file_name, caption) {
  paste0("<figure><img src=\"../plots/", html_escape(file_name), "\" alt=\"",
         html_escape(caption), "\"><figcaption>", html_escape(caption),
         "</figcaption></figure>")
}

build_standalone_report <- function(role_tbl, primary_tbl, extended_tbl, weak_tbl,
                                    readiness_tbl) {
  out <- file.path(AUDIT_REPORT_DIR, "IMRS_Manuscript_Readiness_Report.html")
  primary_row <- primary_tbl %>%
    filter(manuscript_role == "primary_acute_validation") %>%
    slice_head(n = 1)
  extended_row <- extended_tbl %>%
    filter(summary_level == "overall_extended_validation") %>%
    slice_head(n = 1)
  primary_sentence <- if (nrow(primary_row) > 0) {
    paste0(
      "Primary acute external validation includes ", primary_row$n_contrasts,
      " contrasts across ", primary_row$gse_ids,
      "; mean Delta IMRSz = ", sprintf("%.3f", primary_row$mean_delta_imrs_z),
      ", median Delta IMRSz = ", sprintf("%.3f", primary_row$median_delta_imrs_z),
      ", and ", sprintf("%.1f", 100 * primary_row$proportion_positive_delta),
      "% of contrasts are positive."
    )
  } else {
    "No primary acute external validation contrasts were available after cleanup."
  }
  extended_sentence <- if (nrow(extended_row) > 0) {
    paste0(
      "Extended transfer includes ", extended_row$n_contrasts,
      " contrasts across ", extended_row$gse_ids,
      "; mean Delta IMRSz = ", sprintf("%.3f", extended_row$mean_delta_imrs_z),
      ", median Delta IMRSz = ", sprintf("%.3f", extended_row$median_delta_imrs_z),
      ", and ", sprintf("%.1f", 100 * extended_row$proportion_positive_delta),
      "% are positive. These are exploratory, not primary acute validation."
    )
  } else {
    "No extended transfer contrasts were available after cleanup."
  }

  css <- paste0(
    "body{font-family:Arial,Helvetica,sans-serif;margin:36px;line-height:1.45;color:#222;}",
    "h1,h2,h3{color:#16324f;} code{background:#f4f4f4;padding:1px 4px;}",
    ".note{background:#f7fbff;border-left:4px solid #2c7fb8;padding:10px 14px;margin:12px 0;}",
    ".warning{background:#fff8eb;border-left:4px solid #f28e2b;padding:10px 14px;margin:12px 0;}",
    ".table-wrap{max-height:520px;overflow:auto;border:1px solid #ddd;margin:12px 0 24px 0;}",
    "table{border-collapse:collapse;font-size:12px;min-width:100%;}",
    "th,td{border:1px solid #ddd;padding:4px 6px;vertical-align:top;}",
    "th{position:sticky;top:0;background:#f0f3f6;}",
    "figure{margin:16px 0 28px 0;} img{max-width:100%;border:1px solid #ddd;}",
    "figcaption{font-size:13px;color:#555;margin-top:6px;}"
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>IMRS Manuscript Readiness Report</title>",
    "<style>", css, "</style></head><body>",
    "<h1>IMRS Manuscript Readiness Report</h1>",
    "<div class=\"note\"><p>The production model uses five locked acute mouse anchor datasets: GSE39129, GSE167521, GSE264344, GSE279372, and GSE279744. Acute timing alone does not define anchor membership, so acute validation datasets remain validation unless they are part of the selected locked anchor set. Extended and secondary-support datasets are biologically informative but are not primary acute validation evidence.</p><p>Delta IMRSz and directionality are the primary metrics. AUC is secondary because many contrasts have small sample sizes. Weak or near-zero IMRS results are interpreted only after paper/PMID/treatment-context audit.</p></div>",
    "<h2>Manuscript Bottom Line</h2>",
    "<p>", html_escape(primary_sentence), "</p>",
    "<p>", html_escape(extended_sentence), "</p>",
    "<div class=\"warning\"><p>Remaining limitations are explicitly retained as WARNING: small n in many contrasts, no clinical reactogenicity endpoint, mouse-heavy dataset composition, late transfer datasets are not strict acute validation, and some metadata fields require manual verification.</p></div>",
    "<h2>Publication Readiness</h2>",
    table_to_html(readiness_tbl, "Publication readiness summary",
                  "manuscript_publication_readiness_summary.tsv"),
    "<h2>Validation Summaries</h2>",
    table_to_html(primary_tbl, "Manuscript primary validation summary",
                  "manuscript_primary_validation_summary.tsv"),
    table_to_html(extended_tbl, "Extended transfer summary",
                  "manuscript_extended_transfer_summary.tsv"),
    plot_html("manuscript_primary_vs_extended_delta.png", "Primary versus extended Delta IMRSz"),
    plot_html("manuscript_dataset_role_forest.png", "Dataset-role forest plot"),
    "<h2>Cleaned Dataset Roles</h2>",
    table_to_html(role_tbl, "Manuscript dataset role table",
                  "manuscript_dataset_role_table.tsv", max_rows = 80),
    "<h2>Weak Dataset Interpretation</h2>",
    table_to_html(weak_tbl, "Manuscript weak dataset interpretation table",
                  "manuscript_weak_dataset_interpretation_table.tsv"),
    "</body></html>"
  )
  writeLines(html, out, useBytes = TRUE)
  out
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  dataset_path <- file.path(AUDIT_RESULTS_DIR, "dataset_audit_table.tsv")
  step09_path <- file.path(PROJECT_ROOT, "05_score", "transfer", "eval",
                           "step09_split_eval.tsv")
  weak_path <- file.path(AUDIT_RESULTS_DIR, "weak_dataset_paper_context_audit.tsv")

  if (!file.exists(dataset_path)) stop("Missing dataset audit table: ", dataset_path)
  if (!file.exists(step09_path)) stop("Missing Step 09 eval table: ", step09_path)

  dataset_audit <- read_tsv(dataset_path, show_col_types = FALSE, progress = FALSE)
  step09 <- read_tsv(step09_path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      time_h = suppressWarnings(as.numeric(time_h)),
      delta_mean_imrs_z = suppressWarnings(as.numeric(delta_mean_imrs_z)),
      auc_imrs_z = suppressWarnings(as.numeric(auc_imrs_z))
    )

  manual <- manual_metadata()

  role_tbl <- dataset_audit %>%
    mutate(
      dataset_id = gse_id,
      original_dataset_type = dataset_type,
      anchor_calibration_external_original = anchor_calibration_external,
      time_h = suppressWarnings(as.numeric(time_h)),
      locked_anchor_dataset = dataset_id %in% locked_anchor_dataset_ids,
      selected_production_discovery_dataset = locked_anchor_dataset,
      row_is_anchor_phase = case_when(
        is.finite(time_h) & time_h <= 24 ~ TRUE,
        is.finite(time_h) & time_h > 24 ~ FALSE,
        TRUE ~ locked_anchor_dataset
      ),
      production_anchor_row = locked_anchor_dataset & row_is_anchor_phase,
      production_discovery_row = production_anchor_row,
      strict_platform_anchor_row = production_anchor_row,
      additional_acute_discovery_support_row = FALSE,
      manuscript_role = assign_manuscript_role(dataset_id, time_h, split_path),
      manuscript_claim_group = unname(claim_group_map[manuscript_role]),
      final_display_group_v2 = unname(role_display_map[manuscript_role]),
      final_display_group_reason = case_when(
        production_anchor_row ~ "production locked anchor dataset and anchor-phase row",
        dataset_id %in% c("GSE119119", "GSE139529", "GSE279743") ~ "acute validation dataset; acute timing alone does not define discovery-set membership",
        dataset_id == "GSE264344" & is.finite(time_h) & time_h > 24 ~ "same GSE as a locked anchor, but this row is >24 h and therefore extended validation",
        dataset_id %in% c("GSE166655", "GSE178313", "GSE314070") ~ "late or context-shifted validation dataset",
        dataset_id %in% c("GSE262515_cell_line", "GSE262515_tissue") ~ "secondary support / weak-context audit; absent from production Step 06A discovery counts",
        TRUE ~ "not assigned to a main manuscript display group"
      ),
      role_display = final_display_group_v2,
      manuscript_group = final_display_group_v2,
      manuscript_interpretation_group = manuscript_role
    ) %>%
    left_join(manual, by = "dataset_id") %>%
    mutate(
      gse_id_clean = coalesce(gse_id_clean, base_gse_id(dataset_id)),
      organism_clean = coalesce(organism_clean, as.character(organism)),
      role_note = pmap_chr(
        list(dataset_id, time_h, role_note_manual),
        ~ role_note_for_row(..1, ..2, ..3)
      )
    ) %>%
    left_join(
      step09 %>%
        select(gse_id, split_id, delta_mean_imrs_z, auc_imrs_z, pass, fail_reason),
      by = c("dataset_id" = "gse_id", "split_id" = "split_id")
    ) %>%
    select(
      dataset_id, gse_id_clean, split_id,
      original_dataset_type, anchor_calibration_external_original,
      locked_anchor_dataset, selected_production_discovery_dataset,
      row_is_anchor_phase, production_anchor_row,
      production_discovery_row, strict_platform_anchor_row,
      additional_acute_discovery_support_row, final_display_group_v2,
      final_display_group_reason,
      manuscript_role, manuscript_claim_group, role_display,
      manuscript_group, manuscript_interpretation_group,
      organism_original = organism, organism_clean,
      tissue, time_h,
      delivery_platform_original = delivery_platform_type,
      delivery_platform_clean,
      control_definition, n_controls, n_delivery,
      PMID_clean, paper_title_clean, DOI_clean,
      general_biology, role_note,
      delta_mean_imrs_z, auc_imrs_z,
      inclusion_rationale, limitation, split_path,
      everything(), -dataset_type, -anchor_calibration_external, -role_note_manual
    ) %>%
    arrange(match(manuscript_role, names(claim_group_map)), dataset_id, time_h, split_id)

  write_tsv_safe(role_tbl, file.path(AUDIT_RESULTS_DIR, "manuscript_dataset_role_table.tsv"))
  write_csv_safe(role_tbl, file.path(AUDIT_RESULTS_DIR, "manuscript_dataset_role_table.csv"))

  eval_role <- step09 %>%
    left_join(
      role_tbl %>%
        select(dataset_id, split_id, manuscript_role, manuscript_claim_group,
               final_display_group_v2, role_display, manuscript_group,
               manuscript_interpretation_group, organism_clean,
               delivery_platform_clean, general_biology, role_note),
      by = c("gse_id" = "dataset_id", "split_id" = "split_id")
    ) %>%
    mutate(
      manuscript_role = replace_na(manuscript_role, "excluded_or_unclear"),
      manuscript_claim_group = replace_na(manuscript_claim_group, "unclear_not_for_main_claim"),
      final_display_group_v2 = replace_na(final_display_group_v2, "Excluded/unclear"),
      role_display = coalesce(role_display, final_display_group_v2),
      manuscript_group = coalesce(manuscript_group, final_display_group_v2),
      manuscript_interpretation_group = coalesce(manuscript_interpretation_group, manuscript_role)
    )

  role_levels <- names(claim_group_map)
  primary_tbl <- map_dfr(role_levels, function(role) {
    df <- eval_role %>%
      filter(manuscript_role == role, pass %in% TRUE)
    safe_summary(df) %>%
      mutate(
        manuscript_role = role,
        manuscript_claim_group = unname(claim_group_map[role]),
        final_display_group_v2 = unname(role_display_map[role]),
        is_primary_external_validation = role == "primary_acute_validation",
        manuscript_use = case_when(
          role == "anchor" ~ "locked anchor discovery evidence",
          role == "primary_acute_validation" ~ "primary acute validation",
          role == "extended_validation" ~ "extended validation / context-shifted transfer",
          role == "secondary_support" ~ "secondary support / weak-context audit",
          TRUE ~ "not for main claim"
        ),
        .before = 1
      )
  })
  write_tsv_safe(primary_tbl, file.path(AUDIT_RESULTS_DIR,
                                        "manuscript_primary_validation_summary.tsv"))

  extended_overall <- eval_role %>%
    filter(manuscript_role == "extended_validation", pass %in% TRUE) %>%
    safe_summary() %>%
    mutate(summary_level = "overall_extended_validation", dataset_id = "ALL", .before = 1)
  extended_by_dataset <- eval_role %>%
    filter(manuscript_role == "extended_validation", pass %in% TRUE) %>%
    group_by(dataset_id = gse_id) %>%
    group_modify(~ safe_summary(.x)) %>%
    ungroup() %>%
    mutate(summary_level = "by_dataset", .before = 1)
  extended_tbl <- bind_rows(extended_overall, extended_by_dataset) %>%
    select(summary_level, dataset_id, everything())
  write_tsv_safe(extended_tbl, file.path(AUDIT_RESULTS_DIR,
                                         "manuscript_extended_transfer_summary.tsv"))

  if (file.exists(weak_path)) {
    weak_raw <- read_tsv(weak_path, show_col_types = FALSE, progress = FALSE)
    weak_tbl <- weak_raw %>%
      left_join(
        role_tbl %>%
          select(dataset_id, split_id, manuscript_role, manuscript_claim_group,
                 final_display_group_v2, role_display, manuscript_group,
                 manuscript_interpretation_group),
        by = c("dataset_id", "split_id")
      ) %>%
      mutate(
        risk_rank = unname(risk_levels[reviewer_risk_level]),
        risk_rank = replace_na(risk_rank, 0L)
      ) %>%
      select(
        dataset_id, split_id, PMID, paper_title, organism, tissue, time_h,
        treatment_group, control_group, delivery_platform,
        original_IMRS_delta, explanation_category, explainable_status,
        reviewer_risk_level, manuscript_ready_interpretation, notes,
        manuscript_role, manuscript_claim_group
        , final_display_group_v2, role_display, manuscript_group,
        manuscript_interpretation_group
      ) %>%
      arrange(desc(unname(risk_levels[reviewer_risk_level])), dataset_id, time_h)
  } else {
    weak_tbl <- tibble()
  }
  write_tsv_safe(weak_tbl, file.path(AUDIT_RESULTS_DIR,
                                     "manuscript_weak_dataset_interpretation_table.tsv"))

  loao <- read_tsv(file.path(AUDIT_RESULTS_DIR, "leave_one_anchor_out_summary.tsv"),
                   show_col_types = FALSE, progress = FALSE)
  logo <- read_tsv(file.path(AUDIT_RESULTS_DIR, "leave_one_gene_out_summary.tsv"),
                   show_col_types = FALSE, progress = FALSE)
  dominance <- read_tsv(file.path(AUDIT_RESULTS_DIR, "gene_dominance_summary.tsv"),
                        show_col_types = FALSE, progress = FALSE)
  threshold <- read_tsv(file.path(AUDIT_RESULTS_DIR, "threshold_sensitivity_summary.tsv"),
                        show_col_types = FALSE, progress = FALSE)

  anchor_row <- primary_tbl %>% filter(manuscript_role == "anchor") %>% slice_head(n = 1)
  primary_row <- primary_tbl %>% filter(manuscript_role == "primary_acute_validation") %>% slice_head(n = 1)
  extended_row <- extended_tbl %>% filter(summary_level == "overall_extended_validation") %>% slice_head(n = 1)

  logo_preserved <- mean(logo$direction_preserved %in% TRUE, na.rm = TRUE)
  dominance_med <- median(dominance$median_max_contribution_fraction, na.rm = TRUE)
  dominance_max <- max(dominance$max_max_contribution_fraction, na.rm = TRUE)
  threshold_min_prop <- min(threshold$external_proportion_positive_delta, na.rm = TRUE)
  weak_has_warning <- nrow(weak_tbl) > 0 &&
    any(weak_tbl$explainable_status != "explainable" | weak_tbl$reviewer_risk_level != "low")

  readiness_tbl <- tibble::tribble(
    ~evidence_area, ~status, ~manuscript_interpretation, ~reviewer_risk, ~action_needed_before_submission,
    "Anchor derivation evidence",
    ifelse(nrow(anchor_row) > 0 && anchor_row$proportion_positive_delta[[1]] >= 0.90, "PASS", "WARNING"),
    paste0("The five locked anchor datasets are derivation evidence within the production discovery/anchor set; cleaned tables keep them separate from validation. Positive-anchor proportion = ", ifelse(nrow(anchor_row) > 0, sprintf("%.3f", anchor_row$proportion_positive_delta[[1]]), "NA"), "."),
    "low",
    "State that the production locked anchor set is GSE39129, GSE167521, GSE264344, GSE279372, and GSE279744.",
    "Leave-one-anchor-out portability",
    ifelse(nrow(loao) > 0 && all(loao$all_directions_positive %in% TRUE), "PASS", "FAIL"),
    "Held-out anchor analyses retain positive Delta IMRSz in the available anchor-perturbation checks.",
    "low",
    "Report Delta IMRSz as primary; keep AUC secondary.",
    "Leave-one-gene-out non-degeneracy",
    ifelse(is.finite(logo_preserved) && logo_preserved >= 0.95, "PASS", ifelse(is.finite(logo_preserved) && logo_preserved >= 0.85, "WARNING", "FAIL")),
    paste0("Direction is preserved in ", sprintf("%.1f", 100 * logo_preserved), "% of leave-one-gene-out tests."),
    risk_from_status(ifelse(is.finite(logo_preserved) && logo_preserved >= 0.95, "PASS", ifelse(is.finite(logo_preserved) && logo_preserved >= 0.85, "WARNING", "FAIL"))),
    "Mention that high-weight genes were tested for dominance and removal sensitivity.",
    "Gene dominance",
    ifelse(is.finite(dominance_med) && dominance_med < 0.50 && dominance_max < 0.75, "PASS", "WARNING"),
    paste0("Median top-gene contribution fraction across contrast summaries = ", sprintf("%.3f", dominance_med), "; maximum observed = ", sprintf("%.3f", dominance_max), "."),
    ifelse(dominance_max < 0.75, "low", "medium"),
    "Include dominance plot/table as robustness evidence.",
    "Threshold sensitivity",
    ifelse(is.finite(threshold_min_prop) && threshold_min_prop >= 0.75, "PASS", ifelse(is.finite(threshold_min_prop) && threshold_min_prop >= 0.50, "WARNING", "FAIL")),
    paste0("Across threshold grid, minimum external positive-delta proportion = ", sprintf("%.3f", threshold_min_prop), "."),
    risk_from_status(ifelse(is.finite(threshold_min_prop) && threshold_min_prop >= 0.75, "PASS", ifelse(is.finite(threshold_min_prop) && threshold_min_prop >= 0.50, "WARNING", "FAIL"))),
    "Describe parameter grid and report that conclusions are not threshold-fragile.",
    "Primary acute external validation",
    ifelse(nrow(primary_row) > 0 && primary_row$proportion_positive_delta[[1]] >= 0.75, "PASS", ifelse(nrow(primary_row) > 0 && primary_row$proportion_positive_delta[[1]] >= 0.50, "WARNING", "FAIL")),
    paste0("Primary validation is restricted to acute validation contrasts outside the five selected discovery datasets: n = ", ifelse(nrow(primary_row) > 0, primary_row$n_contrasts[[1]], 0), ", mean Delta IMRSz = ", ifelse(nrow(primary_row) > 0, sprintf("%.3f", primary_row$mean_delta_imrs_z[[1]]), "NA"), ", positive proportion = ", ifelse(nrow(primary_row) > 0, sprintf("%.3f", primary_row$proportion_positive_delta[[1]]), "NA"), "."),
    risk_from_status(ifelse(nrow(primary_row) > 0 && primary_row$proportion_positive_delta[[1]] >= 0.75, "PASS", ifelse(nrow(primary_row) > 0 && primary_row$proportion_positive_delta[[1]] >= 0.50, "WARNING", "FAIL"))),
    "Do not pool locked anchors, secondary-support, or extended validation datasets into this primary claim.",
    "Extended transfer interpretation",
    "WARNING",
    paste0("Extended transfer is biologically informative but exploratory: n = ", ifelse(nrow(extended_row) > 0, extended_row$n_contrasts[[1]], 0), ", mean Delta IMRSz = ", ifelse(nrow(extended_row) > 0, sprintf("%.3f", extended_row$mean_delta_imrs_z[[1]]), "NA"), ", positive proportion = ", ifelse(nrow(extended_row) > 0, sprintf("%.3f", extended_row$proportion_positive_delta[[1]]), "NA"), "."),
    "medium",
    "Label late timepoints as exploratory transfer, not primary acute validation.",
    "Weak dataset paper-context audit",
    ifelse(weak_has_warning, "WARNING", "PASS"),
    "Weak or near-zero IMRS results are interpreted only after paper/PMID/treatment-context audit; partially explainable rows remain reviewer-sensitive.",
    ifelse(weak_has_warning, "medium", "low"),
    "Verify paper/treatment metadata for medium-risk rows before manuscript submission.",
    "Remaining limitations",
    "WARNING",
    "Small n in many contrasts, no clinical reactogenicity endpoint, mouse-heavy dataset composition, late transfer datasets are not primary acute validation, and some metadata fields require manual verification.",
    "medium",
    "Keep these limitations explicit in the Discussion and Methods."
  )
  write_tsv_safe(readiness_tbl, file.path(AUDIT_RESULTS_DIR,
                                          "manuscript_publication_readiness_summary.tsv"))

  role_plot_levels <- c(
    "anchor",
    "primary_acute_validation",
    "extended_validation",
    "secondary_support"
  )
  plot_roles <- eval_role %>%
    filter(manuscript_role %in% role_plot_levels,
           is.finite(delta_mean_imrs_z)) %>%
    mutate(
      manuscript_role = factor(manuscript_role,
                               levels = role_plot_levels,
                               labels = unname(role_display_map[role_plot_levels])),
      manuscript_claim_group = factor(manuscript_claim_group,
                                      levels = unname(claim_group_map[role_plot_levels]))
    )
  role_colors <- c(
    "Locked anchor" = "#2F6B9A",
    "Primary acute validation" = "#C84B31",
    "Extended validation" = "#7B5EA7",
    "Secondary support" = "#59A14F"
  )
  p1 <- ggplot(plot_roles, aes(x = manuscript_role, y = delta_mean_imrs_z,
                               fill = manuscript_role)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.58) +
    geom_jitter(aes(color = manuscript_role), width = 0.15, height = 0,
                size = 1.8, alpha = 0.75) +
    scale_fill_manual(values = role_colors, drop = FALSE) +
    scale_color_manual(values = role_colors, drop = FALSE) +
    labs(
      title = "Delivery-associated IMRS shifts by manuscript evidence group",
      x = NULL,
      y = "Delivery-minus-control IMRS z-score",
      fill = "Manuscript role",
      color = "Manuscript role"
    ) +
    theme_audit() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_plot(p1, file.path(AUDIT_PLOTS_DIR,
                          "manuscript_primary_vs_extended_delta.png"),
            width = 8.5, height = 5.5)

  forest_tbl <- plot_roles %>%
    group_by(gse_id, manuscript_role) %>%
    summarise(
      mean_delta_imrs_z = mean(delta_mean_imrs_z, na.rm = TRUE),
      median_delta_imrs_z = median(delta_mean_imrs_z, na.rm = TRUE),
      n_contrasts = n(),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(gse_id, " (n=", n_contrasts, ")"),
      label = factor(label, levels = label[order(mean_delta_imrs_z)])
    )
  p2 <- ggplot(forest_tbl, aes(x = mean_delta_imrs_z, y = label,
                               color = manuscript_role)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 2.5) +
    scale_color_manual(values = role_colors, drop = FALSE) +
    labs(
      title = "Dataset-Level Mean Delta by Manuscript Role",
      x = "Mean delivery-minus-control IMRS z-score",
      y = NULL,
      color = "Manuscript role"
    ) +
    theme_audit() +
    theme(axis.text.x = element_text(angle = 0))
  save_plot(p2, file.path(AUDIT_PLOTS_DIR, "manuscript_dataset_role_forest.png"),
            width = 8.5, height = 6)

  report_path <- build_standalone_report(role_tbl, primary_tbl, extended_tbl,
                                         weak_tbl, readiness_tbl)
  log_msg("INFO", "Wrote manuscript cleanup outputs and report: ", report_path)
  invisible(role_tbl)
}

main()
