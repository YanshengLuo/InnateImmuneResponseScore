project_root <- "D:/IMRS_Project"

all_paths <- list.files(project_root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
info <- file.info(all_paths)

df <- data.frame(
  full_path = rownames(info),
  is_dir = info$isdir,
  size_bytes = info$size,
  stringsAsFactors = FALSE
)

df$rel_path <- sub(
  paste0("^", normalizePath(project_root, winslash = "/"), "/?"),
  "",
  normalizePath(df$full_path, winslash = "/")
)

df$ext <- tools::file_ext(df$rel_path)
df$ext[df$is_dir] <- "DIR"

write.table(
  df[order(df$rel_path), ],
  file = file.path(project_root, "project_inventory_full.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df[!df$is_dir, ][order(-df[!df$is_dir, ]$size_bytes), ],
  file = file.path(project_root, "project_largest_files.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)