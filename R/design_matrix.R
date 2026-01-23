library(dplyr)
library(readr)

make_design <- function(df, dataset_id, tissue_name,
                        project_root = "C:/Users/john/Desktop/IMRS_Project") {
  
  design <- df %>%
    filter(tissue == tissue_name) %>%
    transmute(
      sample_id = srr_id,
      tissue = tissue_name,
      delivery_type = factor(trimws(delivery_type),
                             levels = c("baseline", "Ad26", "Ad5")),
      timepoint_hr = as.numeric(timepoint_hr),
      group = factor(paste0(delivery_type, "_", timepoint_hr)),
      batch = "NA"
    )
  
  outpath <- file.path(project_root, "00_metadata",
                       paste0(dataset_id, "_", tissue_name, "_design.tsv"))
  write_tsv(design, outpath)
}


# generate Blood design
make_design(samples, "GSE264344", "Blood")

# check
read_tsv("C:/Users/john/Desktop/IMRS_Project/00_metadata/GSE264344_Blood_design.tsv", show_col_types = FALSE) %>%
  head()

make_design(samples, "GSE264344", "Muscle")
make_design(samples, "GSE264344", "dLN")
