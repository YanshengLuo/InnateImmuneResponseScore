## ============================================================
## RUN_ALL_METADATA_BUILD.R
## Iterate over all GSE* folders and run METADATA_BUILD.R
## Robust summary binding + notes/warnings surfaced
## ============================================================

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(readr)
  library(tibble)
})

project_root <- "C:/Users/john/Desktop/IMRS_Project"
meta_root <- file.path(project_root, "00_metadata")

metadata_build_script <- file.path(
  project_root,
  "Hypergator_scripts",
  "InnateImmuneResponseScore",
  "R",
  "METADATA_BUILD.R"
)

if (!file.exists(metadata_build_script)) {
  stop("METADATA_BUILD.R not found at:\n", metadata_build_script)
}

# Detect GSE* directories (non-recursive)
all_dirs <- list.dirs(meta_root, full.names = TRUE, recursive = FALSE)
gse_dirs <- all_dirs[str_detect(basename(all_dirs), "^GSE\\d+.*$")]

cat("\nDetected datasets:\n")
print(basename(gse_dirs))

if (length(gse_dirs) == 0) {
  stop("No GSE* directories found under ", meta_root)
}

# Collect rows in a list, then bind_rows at the end
results_list <- list()

for (d in gse_dirs) {
  dataset_id <- basename(d)

  cat("\n============================================================\n")
  cat("Running METADATA_BUILD for:", dataset_id, "\n")
  cat("============================================================\n")

  sra_path <- file.path(meta_root, dataset_id, "SraRunTable.csv")
  if (!file.exists(sra_path)) {
    msg <- "SKIP: Missing SraRunTable.csv"
    cat(msg, "\n")

    results_list[[dataset_id]] <- tibble(
      dataset_id = dataset_id,
      status = "SKIP",
      message = msg,
      warnings = "",
      notes = ""
    )
    next
  }

  # Run metadata build in a clean environment
  run_env <- new.env(parent = globalenv())
  run_env$dataset_id <- dataset_id

  warn_msgs <- character(0)

  res <- tryCatch(
    withCallingHandlers(
      {
        source(metadata_build_script, local = run_env)
        list(status = "OK", message = "")
      },
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(status = "ERROR", message = conditionMessage(e))
    }
  )

  # Pull notes from samples file if present
  notes_msg <- ""
  samples_path <- file.path(meta_root, dataset_id, paste0(dataset_id, "_samples.tsv"))
  if (file.exists(samples_path)) {
    smp <- read_tsv(samples_path, show_col_types = FALSE)
    if ("notes" %in% names(smp)) {
      notes_msg <- smp %>%
        transmute(notes = as.character(notes)) %>%
        filter(!is.na(notes), notes != "") %>%
        distinct() %>%
        pull(notes) %>%
        paste(collapse = "; ")
    }
  }

  results_list[[dataset_id]] <- tibble(
    dataset_id = dataset_id,
    status = res$status,
    message = res$message,
    warnings = paste(unique(warn_msgs), collapse = " | "),
    notes = notes_msg
  )
}

results <- bind_rows(results_list) %>%
  arrange(desc(status == "ERROR"), desc(status == "SKIP"), dataset_id)

cat("\n==================== SUMMARY ====================\n")
print(results, n = Inf)
cat("=================================================\n\n")

# Write summary log
summary_path <- file.path(meta_root, "metadata_build_summary.csv")
write.csv(results, summary_path, row.names = FALSE)
cat("Wrote summary CSV to:\n", summary_path, "\n")
