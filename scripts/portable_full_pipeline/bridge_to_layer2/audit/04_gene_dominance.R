#!/usr/bin/env Rscript

this_script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
  value = TRUE)[1]), winslash = "/", mustWork = FALSE)
source(file.path(dirname(this_script), "00_config.R"))

section_name <- "04_gene_dominance"

mode_chr <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

sample_dominance <- function(comp, weights_df) {
  contrib <- sweep(comp$z, 1, comp$wvec[rownames(comp$z)], "*")
  abs_contrib <- abs(contrib)
  denom <- colSums(abs_contrib, na.rm = TRUE)
  frac <- sweep(abs_contrib, 2, denom, "/")
  frac[, !is.finite(denom) | denom <= 0] <- NA_real_
  max_idx <- apply(frac, 2, function(v) {
    if (all(is.na(v))) return(NA_integer_)
    which.max(v)
  })
  top_gene <- ifelse(is.na(max_idx), NA_character_, rownames(frac)[max_idx])
  symbol_map <- weights_df %>%
    select(gene_id, gene_symbol) %>%
    distinct(gene_id, .keep_all = TRUE)
  tibble(
    sample_id = colnames(frac),
    max_contribution_fraction = as.numeric(frac[cbind(max_idx, seq_along(max_idx))]),
    top_contributor_gene_id = top_gene,
    total_abs_contribution = as.numeric(denom)
  ) %>%
    left_join(symbol_map, by = c("top_contributor_gene_id" = "gene_id")) %>%
    rename(top_contributor_symbol = gene_symbol) %>%
    left_join(comp$design %>% select(sample_id, condition_simple), by = "sample_id")
}

summarize_split_dominance <- function(sample_dom, row) {
  split <- read_split_table(row$split_path)
  merged <- split %>%
    select(any_of(c("sample_id", "condition_simple"))) %>%
    inner_join(sample_dom, by = c("sample_id", "condition_simple"))
  if (nrow(merged) == 0) {
    return(tibble(
      gse_id = row$gse_id,
      split_id = row$split_id,
      split_path = row$split_path,
      pass = FALSE,
      fail_reason = "no_sample_overlap"
    ))
  }
  tibble(
    gse_id = row$gse_id,
    split_id = row$split_id,
    split_path = row$split_path,
    audit_role = assign_audit_role(row$gse_id, row$time_h),
    time_group = assign_time_group(row$time_h),
    tissue = row$tissue,
    time_h = suppressWarnings(as.numeric(row$time_h)),
    n_samples = nrow(merged),
    n_controls = sum(merged$condition_simple == "CONTROL"),
    n_delivery = sum(merged$condition_simple == "DELIVERY"),
    mean_max_contribution_fraction = safe_mean(merged$max_contribution_fraction),
    median_max_contribution_fraction = safe_median(merged$max_contribution_fraction),
    max_max_contribution_fraction = max(merged$max_contribution_fraction, na.rm = TRUE),
    control_median_max_contribution_fraction =
      safe_median(merged$max_contribution_fraction[merged$condition_simple == "CONTROL"]),
    delivery_median_max_contribution_fraction =
      safe_median(merged$max_contribution_fraction[merged$condition_simple == "DELIVERY"]),
    dominant_gene_id_mode = mode_chr(merged$top_contributor_gene_id),
    dominant_gene_symbol_mode = mode_chr(merged$top_contributor_symbol),
    n_unique_top_contributors = n_distinct(merged$top_contributor_gene_id, na.rm = TRUE),
    pass = TRUE,
    fail_reason = NA_character_
  )
}

main <- function() {
  log_msg("INFO", "Starting ", section_name)
  weights <- load_frozen_weights()
  eval_tbl <- load_step09_eval(required = TRUE) %>%
    mutate(pass = as.logical(pass)) %>%
    filter(pass %in% TRUE)
  if (nrow(eval_tbl) == 0) stop("No passing Step 09 contrasts available.")

  datasets <- sort(unique(eval_tbl$gse_id))
  sample_rows <- list()
  summary_tbl <- map_dfr(datasets, function(ds) {
    rows <- eval_tbl %>% filter(gse_id == ds)
    comp <- tryCatch(
      score_dataset_components(ds, weights, enforce_min_coverage = TRUE),
      error = function(e) {
        log_msg("WARN", "Skipping gene dominance dataset ", ds, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(comp)) return(tibble())
    dom <- sample_dominance(comp, weights) %>%
      mutate(gse_id = ds)
    sample_rows[[ds]] <<- dom
    map_dfr(seq_len(nrow(rows)), function(i) {
      tryCatch(
        summarize_split_dominance(dom, rows[i, ]),
        error = function(e) {
          log_msg("WARN", "Dominance split failed for ", rows$split_path[i],
                  ": ", conditionMessage(e))
          tibble(
            gse_id = rows$gse_id[i],
            split_id = rows$split_id[i],
            split_path = rows$split_path[i],
            pass = FALSE,
            fail_reason = conditionMessage(e)
          )
        }
      )
    })
  })

  sample_tbl <- bind_rows(sample_rows)
  out_path <- file.path(AUDIT_RESULTS_DIR, "gene_dominance_summary.tsv")
  write_tsv_safe(summary_tbl, out_path)

  if (nrow(sample_tbl) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No gene dominance results") +
      theme_void()
  } else {
    sample_plot_tbl <- sample_tbl %>%
      left_join(eval_tbl %>% select(gse_id, audit_role, time_group) %>% distinct(),
                by = "gse_id") %>%
      mutate(audit_role = replace_na(audit_role, "external"))
    p <- ggplot(sample_plot_tbl, aes(x = max_contribution_fraction,
                                     fill = audit_role)) +
      geom_histogram(binwidth = 0.025, alpha = 0.72, position = "identity") +
      scale_fill_manual(values = c(anchor = "#4E79A7", calibration = "#59A14F",
                                   external = "#E15759"), drop = FALSE) +
      labs(
        title = "Per-Sample Top Gene Contribution Fraction",
        x = "max abs(w * z) / sum abs(w * z)",
        y = "Samples",
        fill = "Audit role"
      ) +
      theme_audit() +
      theme(axis.text.x = element_text(angle = 0))
  }

  save_plot(p, file.path(AUDIT_PLOTS_DIR, "gene_dominance_distribution.png"),
            width = 7.5, height = 5)

  log_msg("INFO", "Wrote gene dominance rows=", nrow(summary_tbl), " to ", out_path)
  invisible(summary_tbl)
}

main()
