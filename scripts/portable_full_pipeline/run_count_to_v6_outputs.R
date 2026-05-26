#!/usr/bin/env Rscript

this_file <- normalizePath(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)
script_root <- dirname(this_file)
source(file.path(script_root, "lib", "portable_pipeline_utils.R"))

args <- imrs_parse_args()
if (isTRUE(args$help)) {
  cat(
    "Usage: Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R",
    " --config config/full_pipeline_config.yml [--dry-run]",
    " [--mode canonical|all_scored|test] [--start-stage 00] [--stop-stage 07]",
    " [--force] [--skip-manuscript]\n"
  )
  quit(status = 0)
}

ctx <- imrs_load_context(args)
required_locked_anchors <- c("GSE167521", "GSE264344", "GSE279372", "GSE279744", "GSE39129")
if (ctx$mode %in% c("canonical", "all_scored")) {
  missing_config_anchors <- setdiff(required_locked_anchors, ctx$locked_anchors)
  if (length(missing_config_anchors) > 0) {
    stop(
      "Mode ", ctx$mode,
      " requires the fixed locked anchor set in configuration: ",
      paste(required_locked_anchors, collapse = ", "),
      ". Missing from datasets.locked_anchors: ",
      paste(missing_config_anchors, collapse = ", "),
      call. = FALSE
    )
  }
}
imrs_create_output_dirs(ctx)
imrs_require(c("yaml", "readr", "dplyr", "tibble"))
log_msg <- imrs_log_function(ctx, "run_count_to_v6_outputs")
manifest_dir <- file.path(ctx$paths$output_root, "manifests")
label <- imrs_model_label(ctx)
rscript <- if (.Platform$OS.type == "windows") file.path(R.home("bin"), "Rscript.exe") else
  file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) stop("Could not find Rscript executable at ", rscript, ".",
                                call. = FALSE)

ported <- file.path(script_root, "original_ported")
scripts <- list(
  `00` = file.path(script_root, "00_preflight_inputs.R"),
  `01` = c(
    file.path(ported, "01_designs", "02_adding_control_to_design.R"),
    file.path(ported, "01_designs", "03_Final_design_file_autogroup.R")
  ),
  `02` = file.path(ported, "02_deseq2_contrasts", "05_DEseq_contrast_delivery_control.R"),
  `03` = c(
    file.path(ported, "03_modeling", "06A_core_gene_set.R"),
    file.path(ported, "03_modeling", "06B_gene_heterogeneity.R"),
    file.path(ported, "03_modeling", "06C_Power_analysis.R"),
    file.path(ported, "03_modeling", "07_weight_estimation.R"),
    file.path(script_root, "compare_generated_frozen_weights.R")
  ),
  `04` = file.path(ported, "04_scoring", "08_score_samples.R"),
  `05` = file.path(ported, "05_evaluation", "09_calibration_evaluation.R"),
  `06` = file.path(script_root, "package_generated_outputs_for_v6_figures.R"),
  `07` = imrs_normalize(ctx$config$manuscript$runner %||% "run_all_manuscript_outputs_v6.R",
                        ctx$repo_root, FALSE)
)

stages <- tibble::tibble(
  stage_id = names(scripts),
  stage_name = c(
    "preflight_inputs",
    "original_design_and_split_preparation",
    "original_deseq2_contrasts",
    "original_frozen_weight_reconstruction",
    "original_frozen_sample_scoring",
    "original_step09_evaluation",
    "package_generated_outputs_for_v6_figures",
    "optional_v6_manuscript_outputs"
  ),
  script_sequence = vapply(scripts, function(x) {
    paste(vapply(x, imrs_rel_or_path, character(1), root = ctx$repo_root),
          collapse = "; ")
  }, character(1)),
  enabled = c(
    isTRUE(ctx$config$execution$run_preflight %||% TRUE),
    isTRUE(ctx$config$execution$run_design_generation %||% TRUE),
    isTRUE(ctx$config$execution$run_deseq2 %||% TRUE),
    isTRUE(ctx$config$execution$run_weight_reconstruction %||% TRUE),
    isTRUE(ctx$config$execution$run_scoring %||% TRUE),
    isTRUE(ctx$config$execution$run_step09 %||% TRUE),
    isTRUE(ctx$config$execution$package_figure_inputs %||% TRUE),
    isTRUE(ctx$config$execution$run_manuscript_outputs %||% FALSE) &&
      !isTRUE(args$skip_manuscript)
  )
)
valid_stages <- stages$stage_id
if (!args$start_stage %in% valid_stages || !args$stop_stage %in% valid_stages) {
  stop("--start-stage and --stop-stage must be one of: ",
       paste(valid_stages, collapse = ", "), call. = FALSE)
}
if (match(args$start_stage, valid_stages) > match(args$stop_stage, valid_stages)) {
  stop("--start-stage must not come after --stop-stage.", call. = FALSE)
}
stages$selected <- seq_len(nrow(stages)) >= match(args$start_stage, valid_stages) &
  seq_len(nrow(stages)) <= match(args$stop_stage, valid_stages)
stages$status <- ifelse(!stages$selected, "outside_requested_range",
                        ifelse(!stages$enabled, "disabled_by_config", "pending"))
stages$started_at <- NA_character_
stages$finished_at <- NA_character_
stages$exit_status <- NA_integer_
stage_log_names <- c(
  `00` = "00_preflight.log",
  `01` = "01_prepare_designs.log",
  `02` = "02_run_deseq2_contrasts.log",
  `03` = "03_reconstruct_frozen_gene_weights.log",
  `04` = "04_score_samples_with_frozen_weights.log",
  `05` = "05_run_step09_evaluation.log",
  `06` = "06_package_v6_figure_inputs.log",
  `07` = "07_run_manuscript_outputs_v6.log"
)

write_stage_status <- function() {
  readr::write_tsv(stages, file.path(manifest_dir, "full_pipeline_stage_status.tsv"))
  readr::write_tsv(stages, file.path(ctx$paths$output_root, "full_pipeline_stage_status.tsv"))
}

configured <- imrs_configured_datasets(ctx)
input_plan <- dplyr::bind_rows(lapply(configured, function(dataset_id) {
  tibble::tibble(
    dataset_id = dataset_id,
    locked_anchor = dataset_id %in% ctx$locked_anchors,
    count_file = imrs_count_path(ctx, dataset_id),
    count_file_exists = file.exists(imrs_count_path(ctx, dataset_id)),
    metadata_file = imrs_metadata_path(ctx, dataset_id),
    metadata_file_exists = file.exists(imrs_metadata_path(ctx, dataset_id)),
    mode = ctx$mode,
    model_label = label
  )
}))
readr::write_tsv(input_plan, file.path(manifest_dir, "dry_run_or_execution_input_plan.tsv"))
planned_outputs <- c(
  "00_preflight/input_inventory.tsv",
  "01_designs/scoring/<DATASET>/<DATASET>_design.tsv",
  "01_designs/splited/<DATASET>_design/*.tsv",
  "04_de/comparison/<DATASET>/deseq2_contrasts/anchor/*.tsv",
  "05_score/anchors/{core_gene_set,support_by_dataset,gene_heterogeneity,gene_power,gene_weights}.tsv",
  "05_score/transfer/scores/*__imrs_scores.tsv",
  "05_score/transfer/qc/*__{qc_summary,top_contributors}.tsv",
  "05_score/transfer/eval/step09_split_{eval,summary,sample_level}.tsv",
  "06_derived_for_figures/data/derived/figure_inputs/*.tsv",
  "layer2_generated_inputs/ (created by run_step09_to_layer2_inputs.R after all_scored)",
  "manuscript_outputs/ (only when enabled)"
)
run_manifest <- tibble::tibble(
  item = c("run_started_at", "config_file", "project_root", "output_root", "mode",
           "model_label", "locked_anchors", "configured_datasets", "dry_run",
           "force_generated_output_overwrite", "start_stage", "stop_stage",
           "released_derived_overwrite", "planned_outputs"),
  value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), ctx$config_file,
            ctx$project_root, ctx$paths$output_root, ctx$mode, label,
            paste(ctx$locked_anchors, collapse = ";"), paste(configured, collapse = ";"),
            as.character(ctx$dry_run), as.character(args$force), args$start_stage,
            args$stop_stage, "FORBIDDEN", paste(planned_outputs, collapse = "; "))
)
readr::write_tsv(run_manifest, file.path(manifest_dir, "full_pipeline_run_manifest.tsv"))
writeLines(c(
  "IMRS original-script portable count-level pipeline session information",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Mode: ", ctx$mode),
  paste0("Model label: ", label),
  capture.output(utils::sessionInfo())
), file.path(manifest_dir, "session_info.txt"), useBytes = TRUE)
write_stage_status()

child_env <- c(
  IMRS_PORTED_OUTPUT_ROOT = ctx$paths$output_root,
  IMRS_PORTED_COUNTS_ROOT = ctx$paths$counts_root,
  IMRS_PORTED_METADATA_ROOT = ctx$paths$verified_metadata_dir,
  IMRS_PORTED_SCORING_DESIGN_ROOT =
    file.path(ctx$paths$output_root, "01_designs", "scoring"),
  IMRS_PORTED_SPLIT_DESIGN_ROOT =
    file.path(ctx$paths$output_root, "01_designs", "splited")
)

run_child <- function(script, child_args, stage_log, env = child_env) {
  if (!file.exists(script)) stop("Pipeline child script not found: ", script, call. = FALSE)
  if (length(env) > 0) {
    prior <- Sys.getenv(names(env), unset = NA_character_)
    on.exit({
      present <- !is.na(prior)
      if (any(present)) do.call(Sys.setenv, as.list(prior[present]))
      if (any(!present)) Sys.unsetenv(names(prior)[!present])
    }, add = TRUE)
    do.call(Sys.setenv, as.list(env))
  }
  output <- system2(rscript, args = shQuote(c(script, child_args)),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (length(output) > 0) {
    cat(paste(output, collapse = "\n"), "\n")
    cat(paste(output, collapse = "\n"), "\n", file = stage_log, append = TRUE)
  }
  if (status != 0L) {
    stop("Child script failed: ", imrs_rel_or_path(script, ctx$repo_root),
         ". Review ", stage_log, ".", call. = FALSE)
  }
}

verify_exists <- function(paths, stage_id) {
  missing <- paths[!file.exists(paths) & !dir.exists(paths)]
  if (length(missing) > 0) {
    stop("Stage ", stage_id, " completed without required output(s): ",
         paste(missing, collapse = "; "), call. = FALSE)
  }
}

prepare_gse262515_alias_inputs <- function(datasets) {
  if (!"GSE262515" %in% datasets) return(datasets)
  aliases <- c(GSE262515_cell_line = "cell_line", GSE262515_tissue = "tissue")
  scoring_root <- file.path(ctx$paths$output_root, "01_designs", "scoring")
  split_root <- file.path(ctx$paths$output_root, "01_designs", "splited")
  base_design <- file.path(scoring_root, "GSE262515", "GSE262515_design.tsv")
  base_split_dir <- file.path(split_root, "GSE262515_design")
  verify_exists(c(base_design, base_split_dir), "01")
  alias_manifest <- list()
  for (alias in names(aliases)) {
    arm <- aliases[[alias]]
    destination <- file.path(scoring_root, alias, paste0(alias, "_design.tsv"))
    curated_source <- file.path(ctx$paths$verified_metadata_dir, "scoring", alias,
                                paste0(alias, "_design.tsv"))
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(curated_source)) {
      if (!file.copy(curated_source, destination, overwrite = TRUE)) {
        stop("Could not stage curated scoring design for ", alias, ".", call. = FALSE)
      }
      design_source <- curated_source
    } else {
      full_design <- readr::read_tsv(base_design, show_col_types = FALSE, progress = FALSE)
      tissue_col <- names(full_design)[tolower(names(full_design)) %in% c("tissue", "tissue_raw")][1]
      if (is.na(tissue_col)) {
        stop("GSE262515 requires curated alias scoring designs or a tissue column to create ",
             "GSE262515_cell_line and GSE262515_tissue staging inputs.", call. = FALSE)
      }
      subset <- full_design[tolower(as.character(full_design[[tissue_col]])) == arm, ,
                            drop = FALSE]
      if (nrow(subset) == 0) {
        stop("No GSE262515 ", arm, " rows found for alias scoring staging.", call. = FALSE)
      }
      readr::write_tsv(subset, destination)
      design_source <- base_design
    }
    alias_split_dir <- file.path(split_root, paste0(alias, "_design"))
    dir.create(alias_split_dir, recursive = TRUE, showWarnings = FALSE)
    candidates <- list.files(base_split_dir, pattern = "\\.tsv$", full.names = TRUE)
    arm_files <- candidates[grepl(paste0("__T=", arm, "__"), basename(candidates),
                                  ignore.case = TRUE)]
    if (length(arm_files) == 0) {
      stop("No generated GSE262515 split contrast files found for ", alias, ".",
           call. = FALSE)
    }
    ok <- file.copy(arm_files, file.path(alias_split_dir, basename(arm_files)),
                    overwrite = TRUE)
    if (!all(ok)) stop("Could not stage split designs for ", alias, ".", call. = FALSE)
    alias_manifest[[length(alias_manifest) + 1L]] <- tibble::tibble(
      evaluation_dataset_id = alias,
      count_dataset_id = "GSE262515",
      scoring_design_source = design_source,
      staged_scoring_design = destination,
      staged_split_dir = alias_split_dir,
      adapter_note = "GSE262515 manuscript arm alias staging; scoring/evaluation algorithm unchanged"
    )
  }
  readr::write_tsv(dplyr::bind_rows(alias_manifest),
                   file.path(manifest_dir, "gse262515_alias_staging_manifest.tsv"))
  unique(c(setdiff(datasets, "GSE262515"), names(aliases)))
}

stage_step08_inputs <- function(datasets) {
  run_token <- format(Sys.time(), "%Y%m%d_%H%M%S")
  staging_root <- file.path(ctx$paths$output_root, "staging", "04_scoring_inputs", run_token)
  staged_counts_root <- file.path(staging_root, "counts")
  staged_design_root <- file.path(staging_root, "scoring_designs")
  entries <- lapply(datasets, function(dataset_id) {
    count_dataset_id <- if (dataset_id %in% c("GSE262515_cell_line", "GSE262515_tissue")) {
      "GSE262515"
    } else {
      dataset_id
    }
    count_source <- imrs_count_path(ctx, count_dataset_id)
    design_source <- file.path(ctx$paths$output_root, "01_designs", "scoring",
                               dataset_id, paste0(dataset_id, "_design.tsv"))
    count_target <- file.path(staged_counts_root, dataset_id, "featurecounts",
                              "validation", "gene_counts_clean.tsv")
    design_target <- file.path(staged_design_root, dataset_id,
                               paste0(dataset_id, "_design.tsv"))
    verify_exists(c(count_source, design_source), "04")
    dir.create(dirname(count_target), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(design_target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(count_source, count_target, overwrite = TRUE)) {
      stop("Could not stage Step08 count input: ", count_source, call. = FALSE)
    }
    if (!file.copy(design_source, design_target, overwrite = TRUE)) {
      stop("Could not stage Step08 scoring design: ", design_source, call. = FALSE)
    }
    tibble::tibble(
      dataset_id = dataset_id,
      count_dataset_id = count_dataset_id,
      staged_input_type = c("counts", "scoring_design"),
      source_path = c(count_source, design_source),
      staged_path = c(count_target, design_target)
    )
  })
  readr::write_tsv(dplyr::bind_rows(entries),
                   file.path(manifest_dir, "step08_staged_input_manifest.tsv"))
  list(counts_root = staged_counts_root, design_root = staged_design_root)
}

prepare_clean_step08_outputs <- function() {
  output_dirs <- file.path(ctx$paths$output_root, "05_score", "transfer",
                           c("scores", "qc"))
  existing <- unlist(lapply(output_dirs, function(path) {
    if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else character()
  }), use.names = FALSE)
  if (length(existing) == 0) return(invisible(TRUE))
  if (!isTRUE(args$force)) {
    stop(
      "Existing generated Step08 score/QC products are present under output_root. ",
      "Use --force to replace them or select a clean output_root; this prevents stale ",
      "score files from changing Step09 dataset scope.",
      call. = FALSE
    )
  }
  output_prefix <- paste0(imrs_normalize(ctx$paths$output_root, NULL, FALSE), "/")
  for (path in output_dirs) {
    normalized <- imrs_normalize(path, NULL, FALSE)
    if (!startsWith(normalized, output_prefix)) {
      stop("Refusing to clean Step08 generated outputs outside configured output_root: ",
           normalized, call. = FALSE)
    }
    unlink(path, recursive = TRUE, force = TRUE)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  log_msg("Cleared prior generated Step08 score/QC files before scoring to prevent cross-mode carryover.")
  invisible(TRUE)
}

write_all_scored_validation_check <- function(scoring_datasets) {
  if (ctx$mode != "all_scored") return(invisible(NULL))
  score_dir <- file.path(ctx$paths$output_root, "05_score", "transfer", "scores")
  eval_path <- file.path(ctx$paths$output_root, "05_score", "transfer", "eval",
                         "step09_split_eval.tsv")
  score_files <- list.files(score_dir, pattern = "__imrs_scores\\.tsv$", full.names = TRUE)
  eval_tbl <- readr::read_tsv(eval_path, show_col_types = FALSE, progress = FALSE)
  primary <- c("GSE119119", "GSE139529", "GSE279743")
  check <- tibble::tibble(
    dataset_id = scoring_datasets,
    count_input_expected = TRUE,
    score_file_exists = vapply(scoring_datasets, function(id) {
      any(grepl(paste0("^", id, "__"), basename(score_files)))
    }, logical(1)),
    step09_rows_present = vapply(scoring_datasets, function(id) {
      any(eval_tbl$gse_id == id)
    }, logical(1)),
    primary_acute_expected = scoring_datasets %in% primary
  )
  readr::write_tsv(check, file.path(manifest_dir, "all_scored_validation_check.tsv"))
  if (any(!check$score_file_exists)) {
    stop("all_scored validation check failed: missing Step08 score files for ",
         paste(check$dataset_id[!check$score_file_exists], collapse = ", "), call. = FALSE)
  }
  missing_primary <- intersect(primary, scoring_datasets)[
    !vapply(intersect(primary, scoring_datasets),
            function(id) any(eval_tbl$gse_id == id), logical(1))
  ]
  if (length(missing_primary) > 0) {
    stop("all_scored validation check failed: Step09 has no rows for primary acute ",
         "validation dataset(s): ", paste(missing_primary, collapse = ", "), call. = FALSE)
  }
  invisible(check)
}

assert_modeling_de_scope <- function(modeling_datasets) {
  comparison_root <- file.path(ctx$paths$output_root, "04_de", "comparison")
  if (!dir.exists(comparison_root)) return(invisible(TRUE))
  dataset_dirs <- list.dirs(comparison_root, recursive = FALSE, full.names = FALSE)
  dataset_dirs <- dataset_dirs[grepl("^GSE", dataset_dirs)]
  unexpected <- setdiff(dataset_dirs, modeling_datasets)
  unexpected_with_anchor_de <- unexpected[vapply(unexpected, function(dataset_id) {
    length(list.files(file.path(comparison_root, dataset_id, "deseq2_contrasts", "anchor"),
                      pattern = "\\.tsv$", full.names = TRUE)) > 0
  }, logical(1))]
  if (length(unexpected_with_anchor_de) > 0) {
    stop(
      "Frozen-weight reconstruction is limited to the five locked anchors, but output_root ",
      "contains non-anchor DE contrast directories from another run: ",
      paste(unexpected_with_anchor_de, collapse = ", "),
      ". Use a clean output_root or remove only those prior generated DE outputs before rerunning.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

log_msg("Original-script portable count-level ", if (ctx$dry_run) "dry run" else "run",
        " started in ", ctx$mode, " mode.")
cat("\nExpected configured inputs\n")
print(input_plan[, c("dataset_id", "locked_anchor", "count_file_exists",
                     "metadata_file_exists", "count_file", "metadata_file")])
discovered_metadata <- if (dir.exists(ctx$paths$verified_metadata_dir)) {
  list.files(ctx$paths$verified_metadata_dir, pattern = "\\.tsv$", full.names = TRUE)
} else character()
cat("\nDiscovered metadata TSV files\n")
cat(if (length(discovered_metadata) > 0) paste(discovered_metadata, collapse = "\n")
    else "(none found)", "\n")
cat("\nPorted original execution order\n")
cat(paste0(" - ", stages$stage_id, ": ", stages$script_sequence), sep = "\n")
cat("\n\nPlanned generated outputs below ", ctx$paths$output_root, "\n", sep = "")
cat(paste0(" - ", planned_outputs), sep = "\n")
cat("\n")

if (!ctx$dry_run && !args$force) {
  existing_product <- file.path(ctx$paths$output_root, "05_score", "anchors", "gene_weights.tsv")
  if (file.exists(existing_product)) {
    stop("Generated pipeline products already exist under output_root. Use --force to overwrite ",
         "generated results, or select a different output_root. Released data/derived remains protected.",
         call. = FALSE)
  }
}

if (ctx$dry_run) {
  stage_log <- file.path(ctx$paths$output_root, "logs", "00_preflight.log")
  run_child(scripts[["00"]], c("--config", ctx$config_file, "--mode", ctx$mode, "--dry-run"),
            stage_log)
  stages$status[stages$stage_id == "00"] <- "completed_dry_run_preflight"
  stages$exit_status[stages$stage_id == "00"] <- 0L
  stages$status[stages$status == "pending"] <- "planned_dry_run_not_executed"
  write_stage_status()
  log_msg("Dry run completed after preflight. DESeq2, modeling, scoring, Step09, and packaging were not executed.")
  quit(status = 0)
}

active_datasets <- NULL
modeling_datasets <- NULL
scoring_datasets <- NULL

for (i in seq_len(nrow(stages))) {
  if (!stages$selected[[i]] || !stages$enabled[[i]]) next
  stage_id <- stages$stage_id[[i]]
  stage_log <- file.path(ctx$paths$output_root, "logs", stage_log_names[[stage_id]])
  writeLines(character(), stage_log)
  stages$status[[i]] <- "running"
  stages$started_at[[i]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  write_stage_status()
  log_msg("Running stage ", stage_id, ": ", stages$stage_name[[i]], ".")
  status <- tryCatch({
    if (stage_id == "00") {
      run_child(scripts[[stage_id]], c("--config", ctx$config_file, "--mode", ctx$mode), stage_log)
    } else {
      if (is.null(active_datasets)) {
        active_datasets <- imrs_active_datasets(ctx)
        modeling_datasets <- if (ctx$mode %in% c("canonical", "all_scored")) {
          ctx$locked_anchors
        } else {
          active_datasets
        }
        scoring_datasets <- active_datasets
      }
      if (stage_id == "01") {
      for (script in scripts[[stage_id]]) run_child(script, ctx$paths$output_root, stage_log)
      verify_exists(file.path(ctx$paths$output_root, "01_designs", "scoring",
                              scoring_datasets, paste0(scoring_datasets, "_design.tsv")), stage_id)
      verify_exists(file.path(ctx$paths$output_root, "01_designs", "splited",
                              paste0(scoring_datasets, "_design")), stage_id)
      scoring_datasets <- prepare_gse262515_alias_inputs(scoring_datasets)
      } else if (stage_id == "02") {
      for (dataset_id in modeling_datasets) {
        dir.create(file.path(ctx$paths$output_root, "04_de", dataset_id),
                   recursive = TRUE, showWarnings = FALSE)
        run_child(scripts[[stage_id]], c(ctx$paths$output_root, dataset_id), stage_log)
      }
      } else if (stage_id == "03") {
      if (ctx$mode == "test" &&
          length(intersect(required_locked_anchors, modeling_datasets)) < length(required_locked_anchors)) {
        stop(
          "Single-anchor test mode cannot reconstruct frozen IMRS weights using the original workflow. ",
          "Provide all locked anchors for canonical reconstruction, or supply an existing ",
          "frozen_gene_weights.tsv for scoring-only test mode.",
          call. = FALSE
        )
      }
      assert_modeling_de_scope(modeling_datasets)
      for (script in scripts[[stage_id]][1:4]) run_child(script, ctx$paths$output_root, stage_log)
      verify_exists(file.path(ctx$paths$output_root, "05_score", "anchors",
                              c("core_gene_set.tsv", "support_by_dataset.tsv",
                                "gene_heterogeneity.tsv", "gene_power.tsv",
                                "gene_weights.tsv")), stage_id)
      run_child(scripts[[stage_id]][[5]], c("--config", ctx$config_file, "--mode", ctx$mode),
                stage_log)
      } else if (stage_id == "04") {
      if ("GSE262515" %in% scoring_datasets) {
        scoring_datasets <- prepare_gse262515_alias_inputs(scoring_datasets)
      }
      prepare_clean_step08_outputs()
      staged_step08 <- stage_step08_inputs(scoring_datasets)
      scoring_env <- child_env
      scoring_env[["IMRS_PORTED_COUNTS_ROOT"]] <- staged_step08$counts_root
      scoring_env[["IMRS_PORTED_SCORING_DESIGN_ROOT"]] <- staged_step08$design_root
      run_child(scripts[[stage_id]], ctx$paths$output_root, stage_log, env = scoring_env)
      score_files <- list.files(file.path(ctx$paths$output_root, "05_score", "transfer", "scores"),
                                pattern = "__imrs_scores\\.tsv$", full.names = TRUE)
      if (length(score_files) < length(scoring_datasets)) {
        stop("Original Step 08 did not produce one score table for every active dataset.",
             call. = FALSE)
      }
      } else if (stage_id == "05") {
      if ("GSE262515" %in% scoring_datasets) {
        scoring_datasets <- prepare_gse262515_alias_inputs(scoring_datasets)
      }
      run_child(scripts[[stage_id]], ctx$paths$output_root, stage_log)
      verify_exists(file.path(ctx$paths$output_root, "05_score", "transfer", "eval",
                              c("step09_split_eval.tsv", "step09_split_summary.tsv",
                                "step09_split_sample_level.tsv")), stage_id)
      write_all_scored_validation_check(scoring_datasets)
      } else if (stage_id == "06") {
      run_child(scripts[[stage_id]], c("--config", ctx$config_file, "--mode", ctx$mode), stage_log)
      verify_exists(file.path(ctx$paths$generated_derived_root, "figure_inputs",
                              "step09_split_eval.tsv"), stage_id)
      } else if (stage_id == "07") {
      staged_config <- file.path(manifest_dir, "manuscript_generated_inputs_config.yml")
      verify_exists(staged_config, stage_id)
      run_child(scripts[[stage_id]], c("--config", staged_config), stage_log, env = character())
      }
    }
    TRUE
  }, error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n", file = stage_log, append = TRUE)
    log_msg("Stage ", stage_id, " failed: ", conditionMessage(e))
    FALSE
  })
  stages$finished_at[[i]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  stages$exit_status[[i]] <- if (status) 0L else 1L
  stages$status[[i]] <- if (status) "completed" else "failed"
  write_stage_status()
  if (!status) stop("Stage ", stage_id, " failed; review ", stage_log, ".",
                    call. = FALSE)
}
log_msg("Original-script portable count-level pipeline completed. Released data/derived inputs were not overwritten.")
