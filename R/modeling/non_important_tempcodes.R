split_dir <- "D:/IMRS_Project/00_metadata/verified_metadata/splited/GSE139529_design"

files <- list.files(split_dir, pattern="\\.tsv$", full.names=TRUE, recursive=TRUE)

for (f in files) {
  b <- basename(f)
  if (!grepl("__VS=", b)) {
    new_b <- sub("\\.tsv$", "__VS=none.tsv", b)
    file.rename(f, file.path(dirname(f), new_b))
  }
}