options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

imrs_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Install them before running this pipeline.",
      call. = FALSE
    )
  }
}

imrs_script_file <- function() {
  candidates <- character()
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) candidates <- c(candidates, sub("^--file=", "", file_arg[[1L]]))
  frame_paths <- unlist(lapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.null(path) || length(path) == 0L) character() else as.character(path[[1L]])
  }), use.names = FALSE)
  candidates <- c(candidates, frame_paths)
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    candidates <- c(candidates, tryCatch(rstudioapi::getActiveDocumentContext()$path,
                                         error = function(e) ""))
  }
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (length(candidates) > 0L) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

imrs_repo_root <- function() {
  env_root <- Sys.getenv("IMRS_REPOSITORY_ROOT", unset = "")
  starts <- unique(c(env_root, dirname(imrs_script_file()), getwd()))
  starts <- starts[!is.na(starts) & nzchar(starts)]
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    if (file.exists(current) && !dir.exists(current)) current <- dirname(current)
    repeat {
      if (file.exists(file.path(current, "README.md")) &&
          dir.exists(file.path(current, "scripts", "portable_full_pipeline"))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not identify the repository root. Keep the extracted repository structure intact.", call. = FALSE)
}

imrs_is_absolute <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]{2}|^/", path)
}

imrs_normalize <- function(path, base = NULL, must_work = FALSE) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  out <- if (imrs_is_absolute(path) || is.null(base)) path else file.path(base, path)
  normalizePath(out, winslash = "/", mustWork = must_work)
}

imrs_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(config = NULL, mode = NULL, dry_run = FALSE,
              start_stage = "00", stop_stage = "07", force = FALSE,
              skip_manuscript = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--dry-run", "--force", "--skip-manuscript")) {
      key <- sub("^--", "", arg)
      key <- gsub("-", "_", key)
      out[[key]] <- TRUE
      i <- i + 1L
    } else if (arg %in% c("--config", "--mode", "--start-stage", "--stop-stage")) {
      if (i == length(args)) stop("Missing value after ", arg, ".", call. = FALSE)
      key <- sub("^--", "", arg)
      key <- gsub("-", "_", key)
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else if (arg %in% c("-h", "--help")) {
      out$help <- TRUE
      i <- i + 1L
    } else {
      stop("Unrecognized argument: ", arg, call. = FALSE)
    }
  }
  out
}

imrs_resolve_config_file <- function(config_arg, repo_root) {
  config_arg <- config_arg %||% file.path("config", "full_pipeline_config.yml")
  candidates <- unique(c(
    imrs_normalize(config_arg, getwd(), FALSE),
    imrs_normalize(config_arg, repo_root, FALSE)
  ))
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop("Configuration file not found: ", config_arg,
         ". Run from the repository root or provide an existing --config path.",
         call. = FALSE)
  }
  hits[[1]]
}

imrs_load_context <- function(args = imrs_parse_args()) {
  imrs_require(c("yaml", "readr", "dplyr", "tibble", "stringr"))
  repo_root <- imrs_repo_root()
  config_file <- imrs_resolve_config_file(args$config, repo_root)
  config <- yaml::read_yaml(config_file)
  mode <- tolower(as.character(args$mode %||% config$mode %||% "canonical"))
  if (!mode %in% c("canonical", "all_scored", "test")) {
    stop("mode must be one of 'canonical', 'all_scored', or 'test'.",
         call. = FALSE)
  }
  project_root <- imrs_normalize(config$project_root %||% ".", repo_root, FALSE)
  required_path_fields <- c("counts_root", "verified_metadata_dir", "output_root",
                            "released_derived_root", "generated_derived_root")
  missing_fields <- setdiff(required_path_fields, names(config$paths %||% list()))
  if (length(missing_fields) > 0) {
    stop("Config paths section is missing: ", paste(missing_fields, collapse = ", "),
         call. = FALSE)
  }
  paths <- lapply(config$paths, imrs_normalize, base = project_root, must_work = FALSE)
  if (mode == "test") {
    paths$output_root <- file.path(paths$output_root, "test_only")
    paths$generated_derived_root <- file.path(paths$output_root, "06_derived_for_figures",
                                              "data", "derived")
  }
  locked <- unlist(config$datasets$locked_anchors %||% character(), use.names = FALSE)
  validation <- unlist(config$datasets$validation %||% locked, use.names = FALSE)
  optional <- unlist(config$datasets$optional_validation %||% character(), use.names = FALSE)
  if (length(locked) == 0) stop("No locked anchor datasets configured.", call. = FALSE)
  list(
    config = config,
    config_file = config_file,
    repo_root = repo_root,
    project_root = project_root,
    paths = paths,
    mode = mode,
    dry_run = isTRUE(args$dry_run),
    args = args,
    locked_anchors = unique(locked),
    validation = unique(validation),
    optional_validation = unique(optional)
  )
}

imrs_count_path <- function(ctx, dataset_id) {
  pattern <- ctx$config$count_file_pattern %||%
    "{dataset}/featurecounts/validation/gene_counts_clean.tsv"
  relative <- gsub("\\{dataset\\}", dataset_id, pattern)
  imrs_normalize(file.path(ctx$paths$counts_root, relative), NULL, FALSE)
}

imrs_metadata_path <- function(ctx, dataset_id) {
  exact <- file.path(ctx$paths$verified_metadata_dir, paste0(dataset_id, "_design.tsv"))
  if (file.exists(exact)) return(imrs_normalize(exact, NULL, FALSE))
  if (!dir.exists(ctx$paths$verified_metadata_dir)) return(imrs_normalize(exact, NULL, FALSE))
  candidates <- list.files(ctx$paths$verified_metadata_dir, pattern = "\\.tsv$",
                           recursive = TRUE, full.names = TRUE)
  matched <- candidates[tolower(basename(candidates)) ==
                          tolower(paste0(dataset_id, "_design.tsv"))]
  if (length(matched) > 0) return(imrs_normalize(matched[[1]], NULL, FALSE))
  imrs_normalize(exact, NULL, FALSE)
}

imrs_configured_datasets <- function(ctx) {
  unique(c(ctx$validation, ctx$optional_validation, ctx$locked_anchors))
}

imrs_available_datasets <- function(ctx) {
  ids <- imrs_configured_datasets(ctx)
  count_paths <- vapply(ids, function(id) imrs_count_path(ctx, id), character(1))
  metadata_paths <- vapply(ids, function(id) imrs_metadata_path(ctx, id), character(1))
  ids[file.exists(count_paths) & file.exists(metadata_paths)]
}

imrs_active_datasets <- function(ctx, require_inputs = TRUE) {
  available <- imrs_available_datasets(ctx)
  if (ctx$mode %in% c("canonical", "all_scored")) {
    missing_anchors <- setdiff(ctx$locked_anchors, available)
    if (require_inputs && length(missing_anchors) > 0) {
      stop(
        "Mode ", ctx$mode,
        " requires all locked anchor datasets for frozen-weight reconstruction. ",
        "Missing prepared counts or metadata for: ",
        paste(missing_anchors, collapse = ", "),
        ". Provide those inputs or rerun explicitly with --mode test for non-canonical output.",
        call. = FALSE
      )
    }
    if (ctx$mode == "all_scored") {
      return(available)
    }
    main <- ctx$locked_anchors
    if (require_inputs) {
      missing_main <- setdiff(main, available)
      if (length(missing_main) > 0) {
        stop("Canonical locked-anchor inputs are missing for: ",
             paste(missing_main, collapse = ", "), call. = FALSE)
      }
    }
    return(intersect(main, available))
  }
  if (require_inputs && length(available) == 0) {
    stop("Test mode found no dataset with both a prepared count matrix and verified metadata.",
         call. = FALSE)
  }
  available
}

imrs_create_output_dirs <- function(ctx) {
  dirs <- file.path(ctx$paths$output_root, c(
    "00_preflight",
    "01_designs/scoring",
    "01_designs/splited",
    "04_de/normalized",
    "04_de/comparison",
    "05_score/anchors",
    "05_score/transfer/scores",
    "05_score/transfer/qc",
    "05_score/transfer/eval",
    "06_derived_for_figures",
    "manuscript_outputs",
    "logs",
    "manifests"
  ))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  invisible(dirs)
}

imrs_log_function <- function(ctx, stage_id) {
  dir.create(file.path(ctx$paths$output_root, "logs"), recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(ctx$paths$output_root, "logs", paste0(stage_id, ".log"))
  writeLines(character(), log_file)
  function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                   paste(..., collapse = ""))
    message(line)
    cat(line, "\n", file = log_file, append = TRUE)
    invisible(line)
  }
}

imrs_first_named_column <- function(df, candidates) {
  if (length(candidates) == 0) return(NA_character_)
  index <- match(tolower(candidates), tolower(names(df)))
  index <- index[!is.na(index)]
  if (length(index) == 0) NA_character_ else names(df)[index[[1]]]
}

imrs_read_metadata <- function(ctx, dataset_id) {
  path <- imrs_metadata_path(ctx, dataset_id)
  if (!file.exists(path)) stop("Verified metadata file not found for ", dataset_id, ": ", path,
                               call. = FALSE)
  meta <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  sample_candidates <- unlist(ctx$config$metadata$sample_id_candidates %||%
                                c("sample_id", "SampleID", "sample", "run", "geo_accession"),
                              use.names = FALSE)
  sample_col <- imrs_first_named_column(meta, sample_candidates)
  if (is.na(sample_col)) {
    stop("Metadata for ", dataset_id, " has no recognized sample ID column. Found: ",
         paste(names(meta), collapse = ", "), call. = FALSE)
  }
  meta$sample_id <- trimws(as.character(meta[[sample_col]]))
  if (any(!nzchar(meta$sample_id)) || anyDuplicated(meta$sample_id)) {
    stop("Metadata for ", dataset_id, " contains empty or duplicated sample IDs.",
         call. = FALSE)
  }
  meta
}

imrs_label_conditions <- function(meta, ctx, dataset_id) {
  condition_candidates <- unlist(ctx$config$metadata$condition_candidates %||%
                                   c("condition_simple", "condition", "group", "treatment"),
                                 use.names = FALSE)
  source_col <- imrs_first_named_column(meta, condition_candidates)
  if (is.na(source_col)) {
    stop("Metadata for ", dataset_id,
         " has no column that can be mapped to treatment/control labels. Found: ",
         paste(names(meta), collapse = ", "), call. = FALSE)
  }
  value <- trimws(as.character(meta[[source_col]]))
  low <- tolower(value)
  low_norm <- gsub("[^a-z0-9]+", "_", low)
  control_config <- tolower(unlist(ctx$config$metadata$control_labels %||% character(),
                                   use.names = FALSE))
  delivery_config <- tolower(unlist(ctx$config$metadata$delivery_labels %||% character(),
                                    use.names = FALSE))
  condition <- rep(NA_character_, length(value))
  rule <- rep(NA_character_, length(value))

  already <- toupper(value) %in% c("CONTROL", "DELIVERY")
  condition[already] <- toupper(value[already])
  rule[already] <- "existing_condition_simple"

  is_control <- low %in% control_config |
    grepl("control|baseline|vehicle|untreat|mock|saline|pbs|delivery_none|(^|_)none($|_)",
          low_norm)
  condition[is.na(condition) & is_control] <- "CONTROL"
  rule[is.na(rule) & is_control] <- "configured_or_control_pattern"

  is_delivery <- low %in% delivery_config |
    grepl("delivery|treated|treat|vector|nanoparticle|lnp|aav|vacc|inject|virus",
          low_norm)
  condition[is.na(condition) & is_delivery] <- "DELIVERY"
  rule[is.na(rule) & is_delivery] <- "configured_or_delivery_pattern"

  fallback <- isTRUE(ctx$config$metadata$treat_non_control_as_delivery %||% TRUE)
  fill <- is.na(condition) & fallback & !is.na(value) & nzchar(value)
  condition[fill] <- "DELIVERY"
  rule[fill] <- "verified_non_control_fallback"

  meta$condition_simple <- condition
  meta$condition_source <- source_col
  meta$condition_mapping_rule <- rule
  if (any(is.na(condition))) {
    bad <- unique(value[is.na(condition)])
    stop("Could not map CONTROL/DELIVERY labels for ", dataset_id,
         ". Unmapped values in ", source_col, ": ", paste(bad, collapse = ", "),
         ". Add labels to config metadata mapping or curate the metadata.",
         call. = FALSE)
  }
  meta
}

imrs_scoring_design <- function(ctx, dataset_id) {
  meta <- imrs_label_conditions(imrs_read_metadata(ctx, dataset_id), ctx, dataset_id)
  meta$dataset_id <- dataset_id
  front <- c("sample_id", "condition_simple", "dataset_id")
  useful <- c("treatment", "group", "timepoint", "timepoint_hr", "time_h", "tissue",
              "platform", "batch", "role", "manuscript_role", "split_id", "contrast_id",
              "condition_source", "condition_mapping_rule")
  keep <- unique(c(front, intersect(useful, names(meta)), names(meta)))
  meta[, keep, drop = FALSE]
}

imrs_read_counts <- function(path, sample_ids, strip_versions = FALSE) {
  if (!file.exists(path)) stop("Count file not found: ", path, call. = FALSE)
  tab <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  gene_col <- imrs_first_named_column(tab, c("gene_id", "Geneid", "GeneID", "gene"))
  if (is.na(gene_col)) {
    stop("Count file lacks a recognized gene identifier column: ", path,
         ". Expected one of gene_id, Geneid, GeneID, gene.", call. = FALSE)
  }
  sample_cols <- intersect(as.character(sample_ids), names(tab))
  if (length(sample_cols) == 0) {
    stop("No count-matrix sample columns match metadata in: ", path, call. = FALSE)
  }
  numeric_cols <- lapply(tab[, sample_cols, drop = FALSE], function(x) {
    suppressWarnings(as.numeric(x))
  })
  mat <- as.matrix(as.data.frame(numeric_cols, check.names = FALSE))
  rownames(mat) <- as.character(tab[[gene_col]])
  colnames(mat) <- sample_cols
  if (anyNA(mat) || any(!is.finite(mat))) {
    stop("Count file contains nonnumeric or non-finite values in matched sample columns: ",
         path, call. = FALSE)
  }
  if (strip_versions) rownames(mat) <- sub("\\.\\d+$", "", rownames(mat))
  if (anyDuplicated(rownames(mat))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }
  list(matrix = mat, gene_column = gene_col, sample_columns = sample_cols,
       n_gene_rows = nrow(tab))
}

imrs_assert_raw_counts <- function(mat, path) {
  if (any(mat < 0)) {
    stop("Negative values found in count file: ", path,
         ". Provide raw nonnegative featureCounts-style counts.", call. = FALSE)
  }
  if (any(abs(mat - round(mat)) > 1e-8)) {
    stop("Non-integer values found in count file: ", path,
         ". This pipeline requires raw integer counts, not TPM/CPM/log-normalized data.",
         call. = FALSE)
  }
  invisible(TRUE)
}

imrs_model_label <- function(ctx) {
  if (ctx$mode == "canonical") return("canonical_locked_anchor_reconstruction")
  if (ctx$mode == "all_scored") return("locked_anchor_model_all_available_scored")
  "TEST_ONLY_NON_CANONICAL"
}

imrs_rel_or_path <- function(path, root) {
  normalized_path <- imrs_normalize(path, NULL, FALSE)
  normalized_root <- paste0(imrs_normalize(root, NULL, FALSE), "/")
  if (startsWith(normalized_path, normalized_root)) {
    substring(normalized_path, nchar(normalized_root) + 1L)
  } else {
    normalized_path
  }
}
