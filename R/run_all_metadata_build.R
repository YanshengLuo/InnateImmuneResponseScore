## ============================================================
## RUN_ALL_METADATA_BUILD.R
## Iterate over all GSE* folders and run METADATA_BUILD.R
## ============================================================

library(stringr)

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
gse_dirs <- all_dirs[str_detect(basename(all_dirs), "^GSE\\d+$")]

cat("\nDetected datasets:\n")
print(basename(gse_dirs))

if (length(gse_dirs) == 0) {
  stop("No GSE* directories found under ", meta_root)
}

results <- data.frame(
  dataset_id = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

for (d in gse_dirs) {
  dataset_id <- basename(d)

  cat("\n============================================================\n")
  cat("Running METADATA_BUILD for:", dataset_id, "\n")
  cat("============================================================\n")

  sra_path <- file.path(meta_root, dataset_id, "SraRunTable.csv")
  if (!file.exists(sra_path)) {
    msg <- paste("SKIP: Missing SraRunTable.csv")
    cat(msg, "\n")
    results <- rbind(results, data.frame(dataset_id, "SKIP", msg))
    next
  }

  # Run metadata build in a clean environment
  run_env <- new.env(parent = globalenv())
  run_env$dataset_id <- dataset_id

  res <- tryCatch(
    {
      source(metadata_build_script, local = run_env)
      "OK"
    },
    error = function(e) {
      paste("ERROR:", conditionMessage(e))
    }
  )

  if (res == "OK") {
    results <- rbind(results, data.frame(dataset_id, "OK", ""))
  } else {
    results <- rbind(results, data.frame(dataset_id, "ERROR", res))
  }
}

cat("\n==================== SUMMARY ====================\n")
print(results)
cat("=================================================\n\n")

# Write summary log
summary_path <- file.path(meta_root, "metadata_build_summary.csv")
write.csv(results, summary_path, row.names = FALSE)
cat("Wrote summary CSV to:\n", summary_path, "\n")
