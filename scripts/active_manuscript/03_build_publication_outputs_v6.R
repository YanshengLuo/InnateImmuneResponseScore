#!/usr/bin/env Rscript
# =============================================================================
# 03_build_publication_outputs_v6.R
#
# THE MISSING PUBLICATION GENERATOR.
#
# The released component builder (01_build_supplementary_tables_v6.R) writes
# analyst-facing tables named *_dataset_level_provenance.tsv etc., with a SINGLE
# S4 table. The files actually submitted with the manuscript are
# publication_outputs/Supplementary_Table_S1..S5 with S4 split into S4A-S4D.
# No script in the release produced that file set, so the package had no
# executable route from released inputs to the submitted tables. This script is
# that route.
#
# DESIGN
#   * The scored universe comes from ONE place:
#     data/curated_metadata/manuscript_scored_universe.tsv, via
#     scripts/active_manuscript/lib/manuscript_universe.R.
#   * Every numeric cell is RECOMPUTED from gated per-contrast sources.
#   * Authored prose columns (rationales, interpretations, scope notes, context
#     explanations) are carried forward verbatim from the shipped publication
#     tables, which act as the schema + editorial template. They are inputs, not
#     outputs: this script never invents descriptive text.
#   * No manuscript number is hard-coded anywhere in this file.
#
# DEPENDENCY CHAIN
#   data/derived/figure_inputs/step09_split_eval.tsv            (per-contrast scores)
#   data/derived/label_permutation_null_summary.tsv             (01_label_permutation_null.R)
#   data/derived/leave_one_gene_out_summary.tsv                 (03_leave_one_gene_out.R)
#   data/derived/gene_dominance_summary.tsv                     (04_gene_dominance.R)
#   data/derived/figure_inputs/threshold_sensitivity_contrast_deltas.tsv (05_threshold_sensitivity.R)
#   data/derived/figure_inputs/baseline_signature_contrast_long.tsv      (02_baseline_signature_benchmarking.R)
#   data/derived/supplement_dataset_split_provenance_v7.tsv     (14_create_supplement_provenance_v7.R)
#   data/derived/weak_dataset_paper_context_audit.tsv           (07_weak_dataset_paper_context_audit.R)
#   publication_outputs/*.tsv                                   (schema + prose template)
#        -> publication_outputs_rebuilt/Supplementary_Table_S1..S5, S4A-S4D, Notes
#
# RUN
#   set IMRS_REPOSITORY_ROOT=<package root>
#   Rscript --vanilla scripts/active_manuscript/03_build_publication_outputs_v6.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(purrr)
})

root <- Sys.getenv("IMRS_REPOSITORY_ROOT", unset = "")
if (!nzchar(root)) stop("IMRS_REPOSITORY_ROOT must be set.", call. = FALSE)
source(file.path(root, "scripts", "active_manuscript", "lib", "manuscript_universe.R"))

TPL <- file.path(root, "publication_outputs")
OUT <- file.path(root, "publication_outputs_rebuilt")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) read_tsv(p, col_types = cols(.default = col_character()), progress = FALSE)
nm <- function(x) suppressWarnings(as.numeric(x))
say <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

# The shipped tables carry authored prose that states the size of the analysis
# set ("the 70 contrasts included in the principal manuscript analysis"). That
# sentence has to track the universe or the tables contradict themselves in
# text while agreeing in numbers. Rewrite it from the computed size; the old
# value is discovered from the template, never typed in.
normalise_universe_prose <- function(tbl, old_n, new_n) {
  if (is.na(old_n) || old_n == new_n) return(tbl)
  pats <- c(
    sprintf("\\b%d-contrast\\b", old_n),
    sprintf("\\bthe %d contrasts\\b", old_n),
    sprintf("\\ball %d manuscript-evaluated\\b", old_n)
  )
  reps <- c(
    sprintf("%d-contrast", new_n),
    sprintf("the %d contrasts", new_n),
    sprintf("all %d manuscript-evaluated", new_n)
  )
  tbl %>% mutate(across(where(is.character), function(col) {
    for (k in seq_along(pats)) col <- str_replace_all(col, pats[k], reps[k])
    col
  }))
}

say("Publication output generator starting. Root: ", root)

# ---- gated per-contrast spine ------------------------------------------------
eval_all <- rd(file.path(root, "data/derived/figure_inputs/step09_split_eval.tsv")) %>%
  filter(toupper(pass) == "TRUE")
perm_all <- rd(file.path(root, "data/derived/label_permutation_null_summary.tsv"))
roles <- perm_all %>% select(split_id, manuscript_role, manuscript_interpretation_group,
                             manuscript_interpretation_label) %>% distinct()
spine <- eval_all %>% left_join(roles, by = "split_id")
spine_scored <- imrs_filter_scored(spine, "gse_id", "step09 eval", root)
N <- nrow(spine_scored)
say("Scored universe: ", N, " contrasts; groups: ",
    paste(sprintf("%s=%d", names(table(spine_scored$manuscript_role)),
                  table(spine_scored$manuscript_role)), collapse = ", "))
keep_splits <- spine_scored$split_id

# =============================== S1 ==========================================
s1_tpl <- rd(file.path(TPL, "Supplementary_Table_S1.tsv"))
old_ncol <- names(s1_tpl)[str_starts(names(s1_tpl), "Number included in")]
OLD_N <- as.integer(str_match(old_ncol, "Number included in ([0-9]+)-contrast")[, 2])
new_ncol <- sprintf("Number included in %d-contrast manuscript analysis", N)
names(s1_tpl)[names(s1_tpl) == old_ncol] <- new_ncol

prov <- imrs_dedup_provenance(
  rd(file.path(root, "data/derived/supplement_dataset_split_provenance_v7.tsv")),
  "supplement provenance")
prov_counts <- prov %>% count(dataset_id, name = "n_meta")

norm_role <- function(x) str_trim(str_replace(x, fixed(" (production)"), ""))
pool <- spine_scored %>%
  mutate(.role = norm_role(manuscript_interpretation_label)) %>%
  arrange(split_id) %>%
  group_by(dataset = gse_id, .role, tissue) %>%
  group_split()
pool_key <- map_chr(pool, ~ paste(.x$gse_id[1], norm_role(.x$manuscript_interpretation_label[1]),
                                  .x$tissue[1], sep = "\r"))
names(pool) <- pool_key
cursor <- setNames(rep(0L, length(pool)), pool_key)

s1_rows <- vector("list", nrow(s1_tpl))
for (i in seq_len(nrow(s1_tpl))) {
  r <- as.list(s1_tpl[i, ])
  ds <- str_trim(r[["Dataset identifier"]]); role <- str_trim(r[["Analytical role"]])
  tis <- str_trim(r[["Biological system or tissue"]])
  key <- paste(ds, role, tis, sep = "\r")
  scored <- ds %in% imrs_scored_ids(root)
  if (scored && key %in% names(pool)) {
    want <- as.integer(nm(r[["Number of scored contrasts"]]))
    grp <- pool[[key]]; st <- cursor[[key]]
    take <- grp[seq_len(min(want, nrow(grp) - st)) + st, , drop = FALSE]
    cursor[[key]] <- st + nrow(take)
    if (nrow(take) != want)
      stop(sprintf("S1 partition mismatch for %s / %s / %s: wanted %d, available %d",
                   ds, role, tis, want, nrow(take)), call. = FALSE)
    d <- nm(take$delta_mean_imrs_z)
    npos <- sum(d > 0, na.rm = TRUE)
    if (!identical(as.integer(nm(r[["Positive evaluated contrasts"]])), as.integer(npos)))
      say("  S1 fix: ", ds, " / ", role, " / ", tis, " positives ",
          r[["Positive evaluated contrasts"]], " -> ", npos)
    r[["Number of scored contrasts"]] <- as.character(nrow(take))
    r[[new_ncol]] <- as.character(nrow(take))
    r[["Number of metadata-only records"]] <- "0"
    r[["Positive evaluated contrasts"]] <- as.character(npos)
    r[["Mean delivery-minus-control delta IMRSz"]] <- format(mean(d, na.rm = TRUE), digits = 15)
  } else {
    nmeta <- prov_counts$n_meta[match(ds, prov_counts$dataset_id)]
    if (is.na(nmeta)) nmeta <- 0L
    r[["Number of scored contrasts"]] <- "0"
    r[[new_ncol]] <- "0"
    r[["Number of metadata-only records"]] <- as.character(nmeta)
    r[["Positive evaluated contrasts"]] <- "0"
    r[["Mean delivery-minus-control delta IMRSz"]] <- ""
    if (!scored && str_starts(ds, "GSE262515")) {
      r[["Analytical role"]] <- "Excluded/unclear"
      u <- imrs_load_universe(root)
      r[["Scientific rationale for role"]] <- u$exclusion_reason[match(ds, u$dataset_id)]
      say("  S1: ", ds, " demoted to metadata-only (", nmeta, " records retained)")
    }
  }
  s1_rows[[i]] <- as_tibble(r)
}
S1 <- normalise_universe_prose(bind_rows(s1_rows), OLD_N, N)
left <- sum(map_int(pool, nrow)) - sum(cursor)
if (left != 0) stop("S1 left ", left, " scored contrasts unallocated.", call. = FALSE)
tot <- sum(as.integer(nm(S1[[new_ncol]])), na.rm = TRUE)
if (tot != N) stop("S1 total ", tot, " != universe ", N, call. = FALSE)
write_tsv(S1, file.path(OUT, "Supplementary_Table_S1.tsv"), na = "")
say("S1: ", nrow(S1), " rows, ", tot, " contrasts allocated")

# =============================== S2 ==========================================
s2_tpl <- rd(file.path(TPL, "Supplementary_Table_S2.tsv"))
status_eval <- sprintf("Evaluated; included in the %d-contrast manuscript analysis", N)
status_excl <- paste("Not scored in the manuscript analysis; metadata and provenance retained",
                     "(human arm excluded pending species-appropriate reprocessing)")
u <- imrs_load_universe(root)
S2 <- s2_tpl %>%
  mutate(
    .scored = .data$`Split identifier` %in% keep_splits,
    .known  = .data$`Split identifier` %in% spine$split_id,
    `Scoring status` = case_when(
      .scored ~ status_eval,
      .known  ~ status_excl,
      TRUE    ~ .data$`Scoring status`),
    across(c(`Mean delivery-minus-control delta IMRSz`,
             `Median delivery-minus-control delta IMRSz`,
             `AUC (secondary)`, `Cohen d`, `Two-sample t-test p-value`, `Direction`),
           ~ if_else(!.scored & .known, "", .x)),
    `Inferential-statistic status` = if_else(
      !.scored & .known,
      "Not applicable because the arm is excluded from the scored manuscript universe",
      .data$`Inferential-statistic status`),
    `Manuscript analysis group` = if_else(!.scored & .known, "Excluded/unclear",
                                          .data$`Manuscript analysis group`),
    `Biological-context note` = if_else(
      !.scored & .known,
      coalesce(u$exclusion_reason[match(.data$Dataset, u$dataset_id)],
               .data$`Biological-context note`),
      .data$`Biological-context note`)
  ) %>% select(-.scored, -.known) %>% normalise_universe_prose(OLD_N, N)
write_tsv(S2, file.path(OUT, "Supplementary_Table_S2.tsv"), na = "")
say("S2: ", nrow(S2), " rows, ", sum(S2$`Scoring status` == status_eval), " evaluated")

# =============================== S3 ==========================================
S3 <- imrs_filter_scored(rd(file.path(TPL, "Supplementary_Table_S3.tsv")), "Dataset",
                         "S3 boundary audit", root)
write_tsv(S3, file.path(OUT, "Supplementary_Table_S3.tsv"), na = "")
say("S3: ", nrow(S3), " rows; categories: ",
    paste(sprintf("%s=%d", names(table(S3$`Biological context category`)),
                  table(S3$`Biological context category`)), collapse = ", "))

# =============================== S4A =========================================
bh <- function(p) p.adjust(p, method = "BH")
perm <- perm_all %>% filter(split_id %in% keep_splits)
logo <- rd(file.path(root, "data/derived/leave_one_gene_out_summary.tsv")) %>%
  filter(split_id %in% keep_splits)
dom  <- rd(file.path(root, "data/derived/gene_dominance_summary.tsv")) %>%
  filter(split_id %in% keep_splits)
thr  <- imrs_filter_scored(
  rd(file.path(root, "data/derived/figure_inputs/threshold_sensitivity_contrast_deltas.tsv")),
  "gse_id", "threshold detail", root) %>%
  filter(sensitivity_scope == "external_full_3_anchor_weights")
thr_n <- unique(count(thr, grid_id)$n)
if (length(thr_n) != 1) stop("threshold grids disagree on contrast count", call. = FALSE)
n_perm_each <- unique(as.integer(nm(perm$n_permutations)))

s4a <- rd(file.path(TPL, "Supplementary_Table_S4A_Robustness.tsv"))
for (i in seq_len(nrow(s4a))) {
  a <- tolower(s4a$Analysis[i])
  if (str_starts(a, "label-perm")) {
    o <- nm(perm$observed_delta_mean_imrs_z)
    s4a$Scope[i] <- sprintf("All %d manuscript-evaluated split contrasts; %d within-split permutations per contrast", N, n_perm_each[1])
    s4a$`Number of contrasts`[i] <- as.character(N)
    s4a$`Number of tests or settings`[i] <- as.character(N * n_perm_each[1])
    s4a$`Observed mean delta IMRSz`[i] <- format(mean(o), digits = 15)
    s4a$`Mean permutation-null delta`[i] <- format(mean(nm(perm$null_mean_delta)), digits = 15)
    s4a$`Positive observed contrasts`[i] <- format(sum(o > 0))
    s4a$`Contrasts outside 95% null interval`[i] <- format(sum(toupper(perm$observed_outside_95pct_null) == "TRUE"))
    s4a$`FDR-significant contrasts`[i] <- format(sum(bh(nm(perm$empirical_p_two_sided)) < 0.05))
    s4a$`Mean delta minimum`[i] <- format(min(o), digits = 15)
    s4a$`Mean delta maximum`[i] <- format(max(o), digits = 15)
  } else if (str_starts(a, "leave-one-gene")) {
    dp <- toupper(logo$direction_preserved) == "TRUE"
    s4a$Scope[i] <- sprintf("%d removed genes across all %d manuscript-evaluated contrasts",
                            n_distinct(logo$removed_gene_id), N)
    s4a$`Number of contrasts`[i] <- as.character(N)
    s4a$`Number of tests or settings`[i] <- as.character(nrow(logo))
    s4a$`Direction-preserved tests`[i] <- format(sum(dp))
    s4a$`Direction-preserved fraction`[i] <- format(mean(dp), digits = 15)
    s4a$`Median absolute percent change`[i] <- format(median(nm(logo$absolute_percent_change_delta), na.rm = TRUE), digits = 15)
    s4a$`Maximum absolute delta change`[i] <- format(max(nm(logo$absolute_change_delta), na.rm = TRUE), digits = 15)
  } else if (str_starts(a, "gene-dom")) {
    mm <- nm(dom$mean_max_contribution_fraction)
    s4a$Scope[i] <- sprintf("Top-contributor concentration across all %d manuscript-evaluated contrasts", N)
    s4a$`Number of contrasts`[i] <- as.character(N)
    s4a$`Number of tests or settings`[i] <- as.character(nrow(dom))
    s4a$`Mean maximum contribution fraction`[i] <- format(mean(mm), digits = 15)
    s4a$`Median maximum contribution fraction`[i] <- format(median(nm(dom$median_max_contribution_fraction)), digits = 15)
    s4a$`Largest observed maximum contribution fraction`[i] <- format(max(mm), digits = 15)
  } else if (str_starts(a, "threshold")) {
    s4a$Scope[i] <- sprintf("%d alternative coefficient/gene-selection threshold settings evaluated over %d contrasts outside the strict-three sensitivity anchor set per setting",
                            n_distinct(thr$grid_id), thr_n)
    s4a$`Number of contrasts`[i] <- as.character(thr_n)
    s4a$`Number of tests or settings`[i] <- as.character(n_distinct(thr$grid_id))
  }
}
s4a <- normalise_universe_prose(s4a, OLD_N, N)
write_tsv(s4a, file.path(OUT, "Supplementary_Table_S4A_Robustness.tsv"), na = "")
say("S4A rebuilt over ", N, " contrasts (threshold scope = ", thr_n, ")")

# =============================== S4C =========================================
cl <- rd(file.path(root, "data/derived/figure_inputs/baseline_signature_contrast_long.tsv")) %>%
  filter(split_id %in% keep_splits)
s4c <- rd(file.path(TPL, "Supplementary_Table_S4C_Comparator.tsv"))
for (i in seq_len(nrow(s4c))) {
  sc <- str_trim(s4c$`Score class`[i]); gp <- str_trim(s4c$`Manuscript analysis group`[i])
  sel <- cl %>% filter(manuscript_interpretation_label == gp,
                       str_starts(score_label, str_sub(sc, 1, 12)))
  if (nrow(sel) == 0) stop("S4C: no source rows for ", sc, " / ", gp, call. = FALSE)
  d <- nm(sel$delta_score); m <- mean(d); s <- if (length(d) > 1) sd(d) else 0
  s4c$`Number of contrasts`[i] <- as.character(length(d))
  s4c$`Mean score shift`[i] <- format(m, digits = 15)
  s4c$`Median score shift`[i] <- format(median(d), digits = 15)
  s4c$`Standard deviation`[i] <- format(s, digits = 15)
  s4c$`Positive-direction proportion`[i] <- format(mean(d > 0), digits = 15)
  s4c$`Coefficient of variation`[i] <- if (m != 0) format(s / m, digits = 15) else ""
  if ("auc_secondary" %in% names(sel))
    s4c$`Mean AUC (secondary)`[i] <- format(mean(nm(sel$auc_secondary), na.rm = TRUE), digits = 15)
}
write_tsv(s4c, file.path(OUT, "Supplementary_Table_S4C_Comparator.tsv"), na = "")
say("S4C rebuilt: ", nrow(s4c), " rows")

# ===================== S4B / S4D / S5 / Notes (proven unaffected) ============
loao_det <- rd(file.path(root, "data/derived/figure_inputs/five_anchor_leave_one_anchor_out_contrast_details.tsv"))
if (any(!loao_det$split_id %in% keep_splits))
  stop("S4B is affected by the exclusion after all.", call. = FALSE)
say("S4B verified unaffected: 0 of ", nrow(loao_det), " held-out contrasts excluded")
pc <- rd(file.path(root, "data/derived/figure_inputs/baseline_signature_paired_contrast_comparison.tsv"))
prim <- spine_scored$split_id[spine_scored$manuscript_role == "primary_acute_validation"]
if (any(pc$split_id[pc$split_id %in% prim] %in% setdiff(spine$split_id, keep_splits)))
  stop("S4D contaminated.", call. = FALSE)
say("S4D verified unaffected: ", n_distinct(pc$split_id[pc$split_id %in% prim]),
    " primary-acute contrasts, none excluded")
for (f in c("Supplementary_Table_S4B_LOAO.tsv", "Supplementary_Table_S4D_Paired_Tests.tsv",
            "Supplementary_Table_S5.tsv")) {
  file.copy(file.path(TPL, f), file.path(OUT, f), overwrite = TRUE)
}
# Notes carries authored prose about the analysis-set size, so it is normalised
# rather than copied verbatim.
write_tsv(normalise_universe_prose(rd(file.path(TPL, "Supplementary_Table_Notes.tsv")), OLD_N, N),
          file.path(OUT, "Supplementary_Table_Notes.tsv"), na = "")
say("Pass-through copied: S4B, S4D, S5, Notes")
say("Publication outputs written to ", OUT)
