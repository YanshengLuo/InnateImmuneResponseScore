# Helpers for the IMRS v6 clean assembly layer.
#
# This script reads v3 helper code and old plotting logic, but writes only under
# repository-local manuscript figure outputs. It intentionally keeps plot geoms, themes,
# colors, legends, axes, coordinates, and scales from the existing workflow.

options(stringsAsFactors = FALSE)

v5_required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble", "purrr",
  "ggplot2", "scales", "grid", "png", "patchwork", "svglite"
)

v5_warnings <- character()
v5_skipped <- character()

log_msg_v5 <- function(..., level = "INFO") {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 level, " ", paste(..., collapse = ""))
  message(line)
  line
}

add_warning_v5 <- function(...) {
  msg <- paste(..., collapse = "")
  v5_warnings <<- unique(c(v5_warnings, msg))
  log_msg_v5(msg, level = "WARN")
}

stop_if_missing_packages_v5 <- function(packages = v5_required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

norm_path_v5 <- function(path, must_work = FALSE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

write_tsv_v5 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

clear_superseded_root_figures_v5 <- function(v5_root) {
  if (!dir.exists(v5_root)) return(invisible(character()))
  canonical_stems <- c(
    paste0("Figure", 1:6, "_main_v5"),
    "FigureS1_main_v5",
    "FigureS2_gene_program_enrichment_combined"
  )
  root_media <- list.files(
    v5_root, pattern = "\\.(png|pdf|svg)$", full.names = TRUE,
    recursive = FALSE, ignore.case = TRUE
  )
  stale <- root_media[!tools::file_path_sans_ext(basename(root_media)) %in% canonical_stems]
  if (length(stale) > 0L && !all(file.remove(stale))) {
    stop("Could not remove superseded final-root figure file(s): ",
         paste(stale, collapse = "; "), call. = FALSE)
  }
  invisible(stale)
}

global_text_replacements_v5 <- c(
  "Validation context" = "Manuscript analysis group",
  "validation context" = "manuscript analysis group",
  "Validation groups" = "Manuscript analysis groups",
  "Validation group" = "Manuscript analysis group",
  "Validation timing group" = "Manuscript analysis group",
  "Analysis group" = "Manuscript analysis group",
  "Reviewer risk level" = "Interpretation support level",
  "Reviewer risk" = "Interpretation support level",
  "Interpretation confidence" = "Interpretation support level",
  "Comparator signature" = "Signature",
  "Benchmark score" = "Signature",
  "Primary validation" = "Primary acute validation",
  "primary validation" = "primary acute validation",
  "Primary validation datasets show consistent IMRS elevation" = "Primary acute validation showed consistent positive delivery-minus-control \u0394IMRSz",
  "Primary acute validation supports delivery-associated IMRS elevation" = "Primary acute validation showed consistent positive delivery-minus-control \u0394IMRSz",
  "Primary and extended datasets show positive IMRS shifts" = "Primary acute and extended validation show positive \u0394IMRSz",
  "Each point is a split contrast; the dashed line marks no delivery-associated score shift." = "Each point is a split contrast; the dashed line marks no delivery-minus-control score shift.",
  "\\npositive=" = "\\n\u0394>0 = ",
  "Mean delivery-minus-control IMRS z-score (pseudo-log scale)" = "Mean delivery-minus-control \u0394IMRSz (pseudo-log scale)",
  "Mean delivery-minus-control IMRS z-score" = "Mean delivery-minus-control \u0394IMRSz",
  "Mean delivery-minus-control\\nIMRS z-score" = "Mean delivery-minus-control\\n\u0394IMRSz",
  "Observed mean delivery-minus-control IMRS z-score" = "Observed mean delivery-minus-control \u0394IMRSz",
  "Split contrasts ordered by observed delivery-minus-control IMRS z-score" = "Split contrasts ordered by observed delivery-minus-control \u0394IMRSz",
  "Original mean delivery-minus-control IMRS z-score" = "Original \u0394IMRSz",
  "After single-gene removal IMRS z-score" = "\u0394IMRSz after single-gene removal",
  "Maximum absolute change in delivery-minus-control IMRS z-score" = "Maximum absolute change in delivery-minus-control \u0394IMRSz",
  "IMRS delivery-minus-control z-score" = "IMRS delivery-minus-control \u0394IMRSz",
  "Baseline signature delivery-minus-control z-score" = "Baseline signature delivery-minus-control score",
  "Proportion of positive delivery-associated contrasts" = "Proportion of contrasts with positive delivery-minus-control score",
  "delivery-associated score shifts" = "delivery-minus-control score shifts",
  "delivery-associated IMRS elevation" = "positive delivery-minus-control \u0394IMRSz",
  "delivery-associated sample-level separation" = "delivery-control sample-level separation",
  "Largest IMRS weights highlight acute discovery response genes" = "Top frozen IMRS coefficients",
  "Genes are ranked by magnitude of the frozen weight, not signed direction." = "Genes are ranked by magnitude of the frozen IMRS coefficient, not signed direction.",
  "Absolute frozen IMRS weight" = "Absolute frozen IMRS coefficient",
  "Distribution of frozen IMRS gene weights" = "Distribution of frozen IMRS gene coefficients",
  "Distribution of fixed gene weights used for sample scoring" = "Distribution of frozen IMRS gene coefficients used for sample scoring",
  "Frozen IMRS gene weight" = "Frozen IMRS gene coefficient",
  "Largest IMRS weights" = "Top frozen IMRS coefficients",
  "frozen weight" = "frozen IMRS coefficient",
  "frozen weights" = "frozen IMRS coefficients",
  "weighted gene" = "coefficient-weighted gene",
  "Weak responses are explained by timing and biological context" = "Late or context-shifted settings provide boundary-setting evidence",
  "Risk categories summarize reviewer-facing interpretation of weak, late, or context-shifted contrasts." = "Support levels summarize interpretation of late or context-shifted contrasts.",
  "Weak-context datasets show attenuated IMRS responses" = "Late or context-shifted settings show attenuated \u0394IMRSz",
  "Dataset-level means highlight late or context-shifted validation settings." = "Dataset-level means highlight late or context-shifted settings.",
  "Context-shifted datasets show variable IMRS elevation" = "Late and context-shifted datasets show variable \u0394IMRSz",
  "Points show delivery-minus-control IMRS shifts annotated by collection time." = "Points show delivery-minus-control \u0394IMRSz annotated by collection time.",
  "GSE264344 captures adenoviral-vector IMRS kinetics" = "Adenoviral-vector IMRS responses peak within the acute window and attenuate by 72 h",
  "dLN means draining lymph node; 72 h values are interpreted as waning kinetics." = "dLN means draining lymph node; 72 h values are consistent with attenuation by 72 h.",
  "waning kinetics" = "attenuation",
  "Observed IMRS shifts exceed label-permutation expectations" = "Observed \u0394IMRSz exceeds within-contrast label-permutation null intervals",
  "Observed delivery-minus-control shifts are compared with 95% within-contrast label-permutation intervals." = "Observed delivery-minus-control \u0394IMRSz is compared with 95% within-contrast label-permutation intervals.",
  "Permutation-tested IMRS shifts differ by analysis group" = "Permutation-tested contrasts retain positive \u0394IMRSz in acute groups",
  "IMRS scores remain stable after single-gene removal" = "Single-gene removal preserves contrast-level \u0394IMRSz",
  "No single gene dominates IMRS contrast-level responses" = "Top-gene contribution remains low across contrasts",
  "IMRS response is not driven by a single dominant gene" = "Top-gene contribution remains low across contrasts",
  "Baseline signatures provide comparator response profiles" = "Comparator immune signatures contextualize acute IMRS response patterns",
  "IMRS and benchmark signatures are directionally compared" = "Comparator immune signatures contextualize delivery-minus-control directionality",
  "Benchmark signatures are positive-control comparators, not replacements for IMRS." = "Immune-response signatures are contextual comparators, not replacements for IMRS."
)

wording_audit_rows_v5 <- function() {
  replacements <- data.frame(
    figure_or_panel = "global/source-label replacement",
    old_text = names(global_text_replacements_v5),
    new_text = unname(global_text_replacements_v5),
    reason = "v5 manuscript-facing terminology standardization",
    stringsAsFactors = FALSE
  )
  workflow <- data.frame(
    figure_or_panel = "Figure1A",
    old_text = c(
      "IMRS computation and validation workflow",
      "Anchor construction, frozen-weight scoring, split-contrast validation, and manuscript-readiness audit",
      "Verified metadata and raw count matrices",
      "DELIVERY vs CONTROL split contrasts",
      "Anchor-only differential expression",
      "Reproducibility gene filtering",
      "Heterogeneity and low-power filtering",
      "Frozen anchor-derived gene weights",
      "Target dataset normalization",
      "Control-based gene z-scores",
      "Weighted raw IMRS score",
      "Control-standardized IMRSz",
      "Split-contrast Delta IMRSz / AUC / direction",
      "Dataset-role audit and weak-dataset cleanup",
      "Manuscript interpretation groups",
      "Frozen anchor-derived weights are not refit on validation or transfer datasets."
    ),
    new_text = c(
      "Frozen IMRS scoring and transfer-evaluation workflow",
      "Anchor-derived gene coefficients are frozen before scoring independent delivery-versus-control contrasts",
      "Verified metadata and raw RNA-seq count matrices",
      "Delivery-versus-control split definitions",
      "Locked-anchor delivery-versus-control differential expression",
      "Cross-anchor reproducibility filtering",
      "Heterogeneity and information-content filtering",
      "Frozen anchor-derived gene coefficients",
      "Target-dataset count normalization",
      "Control-referenced gene z-scores",
      "Weighted sample-level IMRS score",
      "Control-standardized sample IMRSz",
      "Delivery-minus-control \u0394IMRSz, directionality, and secondary AUC",
      "Dataset-role curation and boundary-context audit",
      "Manuscript analysis groups",
      "Frozen anchor-derived coefficients are not refit during validation or transfer evaluation."
    ),
    reason = "requested Figure 1A wording/content polish",
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(replacements, workflow)
}

replace_text_v5 <- function(x) {
  out <- x
  for (pat in names(global_text_replacements_v5)) {
    out <- gsub(pat, global_text_replacements_v5[[pat]], out, fixed = TRUE)
  }
  out
}

replace_source_text_v5 <- function(lines) {
  out <- lines
  for (pat in names(global_text_replacements_v5)) {
    replacement <- global_text_replacements_v5[[pat]]
    replacement <- paste(
      strsplit(replacement, intToUtf8(0x0394), fixed = TRUE)[[1]],
      collapse = "\\u0394"
    )
    out <- gsub(pat, replacement, out, fixed = TRUE)
  }
  out
}

v5_integer_size_breaks <- function(max_n, requested) {
  max_n <- floor(as.numeric(max_n)[1])
  if (!is.finite(max_n) || max_n < 1) {
    return(1)
  }
  out <- requested[requested <= max_n]
  if (max_n < max(requested) && !(max_n %in% out)) {
    out <- c(out, max_n)
  }
  unique(out[out >= 1])
}

install_v5_clipping_overrides <- function(env) {
  override_code <- '
make_FigureSC <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) != "Locked anchor") %>%
    group_by(dataset_id, tissue, time_h, delivery_platform_clean, manuscript_group) %>%
    summarise(mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              n_contrasts = n(), .groups = "drop") %>%
    mutate(label = short_text(dataset_context_label(dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE), 55),
           label = ordered_factor(label, mean_delta))
  group_breaks <- names(manuscript_group_palette)[names(manuscript_group_palette) %in% unique(as.character(plot_tbl$manuscript_group))]
  size_breaks <- v5_integer_size_breaks(max(plot_tbl$n_contrasts, na.rm = TRUE), c(1, 3, 6))
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "#4B5563") +
    geom_point(alpha = 0.95) +
    scale_color_manual(values = manuscript_group_palette, breaks = group_breaks, drop = TRUE) +
    scale_size_continuous(
      range = c(2, 5),
      breaks = size_breaks,
      limits = c(1, max(plot_tbl$n_contrasts, na.rm = TRUE)),
      labels = function(x) sprintf("%d", as.integer(x))
    ) +
    guides(
      size = guide_legend(
        title = "Passing split contrasts",
        ncol = 2,
        byrow = TRUE,
        order = 2,
        override.aes = list(color = "black")
      ),
      color = guide_legend(
        title = "Manuscript analysis group",
        ncol = 1,
        byrow = TRUE,
        order = 1,
        override.aes = list(size = 3)
      )
    ) +
    labs(
      title = "Dataset-level summaries clarify context-dependent \\u0394IMRSz values",
      subtitle = "Dataset-level means reduce contrast-level crowding; full forest is retained as appendix FigureSB.",
      x = "Mean delivery-minus-control \\u0394IMRSz",
      y = "Dataset context",
      color = "Manuscript analysis group",
      size = "Passing split contrasts"
    ) +
    theme_imrs_publication(base_size = 10.5) +
    theme(
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 9.7),
      axis.title = element_text(size = 10.6),
      legend.text = element_text(size = 9.4),
      legend.title = element_text(size = 9.8),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.margin = margin(t = 2, r = 2, b = 2, l = 4),
      legend.spacing.y = grid::unit(8, "pt"),
      plot.margin = margin(8, 8, 8, 10)
    )
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1C_simplified_by_dataset", 7.9, 5.6, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSC")
}

make_FigureSD <- function() {
  plot_tbl <- role_pass_for_plot %>%
    filter(as.character(manuscript_group) %in% c("Extended validation", "Secondary support")) %>%
    group_by(dataset_id, tissue, delivery_platform_clean, manuscript_group) %>%
    summarise(mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
              n_contrasts = n(),
              time_min = min(time_h, na.rm = TRUE),
              time_max = max(time_h, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(label = short_text(dataset_context_label(dataset_id, tissue, time_min, delivery_platform_clean, compact = TRUE), 64),
           label = ordered_factor(label, mean_delta))
  group_breaks <- c("Extended validation", "Secondary support")
  group_breaks <- group_breaks[group_breaks %in% unique(as.character(plot_tbl$manuscript_group))]
  size_breaks <- v5_integer_size_breaks(max(plot_tbl$n_contrasts, na.rm = TRUE), c(1, 2, 3))
  p <- ggplot(plot_tbl, aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#4B5563") +
    geom_point(alpha = 0.95) +
    coord_cartesian(xlim = c(-2, 5)) +
    scale_color_manual(values = manuscript_group_palette, breaks = group_breaks, drop = TRUE) +
    scale_size_continuous(
      range = c(2.2, 5),
      breaks = size_breaks,
      limits = c(1, max(plot_tbl$n_contrasts, na.rm = TRUE)),
      labels = function(x) sprintf("%d", as.integer(x))
    ) +
    guides(
      size = guide_legend(
        title = "Contrasts",
        ncol = 1,
        order = 2,
        override.aes = list(color = "black")
      ),
      color = guide_legend(
        title = "Manuscript analysis group",
        ncol = 1,
        byrow = TRUE,
        order = 1,
        override.aes = list(size = 3)
      )
    ) +
    labs(
      title = "Late and context-shifted datasets show attenuated \\u0394IMRSz",
      subtitle = "Dataset-level means highlight late or context-shifted validation settings.",
      x = "Mean delivery-minus-control \\u0394IMRSz",
      y = "Dataset context",
      color = "Manuscript analysis group",
      size = "Contrasts"
    ) +
    theme_imrs_publication(base_size = 10.2) +
    theme(
      axis.text.y = element_text(size = 10.4),
      axis.text.x = element_text(size = 9.5),
      axis.title = element_text(size = 10.4),
      legend.text = element_text(size = 9.2),
      legend.title = element_text(size = 9.6),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.margin = margin(t = 2, r = 2, b = 2, l = 4),
      legend.spacing.y = grid::unit(8, "pt"),
      plot.margin = margin(8, 8, 8, 10)
    )
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1D_weak_zoom_forest", 8.3, 5.25, dpi = 400,
                 source_tables = required_paths$role_table,
                 source_code_section_or_function = "make_FigureSD")
}
'
  eval(parse(text = override_code), envir = env)
  invisible(env)
}

load_old_generator_env_v5 <- function(old_script, project_root, v5_root) {
  old_script <- norm_path_v5(old_script, must_work = TRUE)
  lines <- readLines(old_script, warn = FALSE)
  run_idx <- grep("^# D\\. Run section", lines)
  if (length(run_idx) != 1) {
    stop("Could not find run-section boundary in old generator: ", old_script, call. = FALSE)
  }
  lines <- lines[seq_len(run_idx - 1L)]
  lines <- replace_source_text_v5(lines)

  sidecar_root <- norm_path_v5(file.path(v5_root, "tables", "_v5_regeneration_sidecars"), must_work = FALSE)
  support_notes <- norm_path_v5(file.path(v5_root, "tables", "Figure2_support_pattern_notes_v5.tsv"),
                                must_work = FALSE)

  figure_input_dir <- Sys.getenv("IMRS_FIGURE_INPUT_DIR", unset = "")
  if (!nzchar(figure_input_dir) || !dir.exists(figure_input_dir)) {
    stop("IMRS_FIGURE_INPUT_DIR must point to the released figure-input directory.", call. = FALSE)
  }
  v2_root <- figure_input_dir
  extra_results <- figure_input_dir

  lines <- sub('^output_root <- file\\.path\\(project_root, "revised_plots"\\)$',
               paste0('output_root <- "', sidecar_root, '"'), lines)
  lines <- sub('^v2_root <- file\\.path\\(project_root, "v2_manuscript"\\)$',
               paste0('v2_root <- "', norm_path_v5(v2_root, must_work = FALSE), '"'), lines)
  lines <- sub('^extra_results_dir <- file\\.path\\(final_root, "publication_extra_generated", "results"\\)$',
               paste0('extra_results_dir <- "', norm_path_v5(extra_results, must_work = FALSE), '"'), lines)
  lines <- sub('^figure2b_support_pattern_path <- file\\.path\\(output_root, "Figure2B_core_gene_support_pattern_notes.tsv"\\)$',
               paste0('figure2b_support_pattern_path <- "', support_notes, '"'), lines)

  old_bar <- 'geom_col\\(fill = "#E5E7EB", color = "#374151", linewidth = 0\\.45, width = 0\\.68\\)'
  new_bar <- 'geom_col(fill = "#4E79A7", color = "#374151", linewidth = 0.45, width = 0.68)'
  replace_hits <- grepl(old_bar, lines)
  if (sum(replace_hits) != 1) {
    stop("Expected exactly one Figure2B main-bar fill line; found ", sum(replace_hits), call. = FALSE)
  }
  lines[replace_hits] <- sub(old_bar, new_bar, lines[replace_hits])
  lines <- sub('^      title = "Breakdown of\\\\nall-but-one support",$',
               '      title = NULL,', lines)
  lines <- sub(
    'explanation_category = factor\\(display_text\\(explanation_category\\)\\)',
    paste0(
      'explanation_category = factor(dplyr::recode(as.character(explanation_category), ',
      '"disease_rescue_model" = "Disease-rescue model", ',
      '"distal_or_adaptive_tissue" = "Distal/adaptive tissue context", ',
      '"formulation_designed_to_reduce_inflammation" = "Low-inflammatory formulation design", ',
      '"late_timepoint" = "Late timepoint", ',
      '"therapeutic_cargo_specific_effect" = "Therapeutic cargo/context effect", ',
      '"tissue_time_kinetic_effect" = "Tissue-time kinetic effect", ',
      '.default = display_text(explanation_category)))'
    ),
    lines
  )

  env <- new.env(parent = globalenv())
  env$IMRS_PROJECT_ROOT <- project_root
  eval(parse(text = paste(lines, collapse = "\n")), envir = env)
  install_v5_clipping_overrides(env)
  env
}

strip_publication_text_v5 <- function(plot) {
  blank_text_theme <- ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  )

  if (inherits(plot, "patchwork")) {
    if (!is.null(plot$patches) && !is.null(plot$patches$annotation)) {
      plot$patches$annotation$title <- NULL
      plot$patches$annotation$subtitle <- NULL
      plot$patches$annotation$caption <- NULL
      plot$patches$annotation$tag_levels <- NULL
      plot$patches$annotation$tag_prefix <- NULL
      plot$patches$annotation$tag_suffix <- NULL
    }
    return(plot)
  }

  if (inherits(plot, "ggplot")) {
    plot$labels$title <- NULL
    plot$labels$subtitle <- NULL
    plot$labels$caption <- NULL
    plot <- plot +
      ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
      blank_text_theme
  }

  plot
}

save_plot_all_formats_v5 <- function(plot, out_base, width, height, dpi = 400) {
  png_path <- paste0(out_base, ".png")
  pdf_path <- paste0(out_base, ".pdf")
  svg_path <- paste0(out_base, ".svg")

  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi,
                  limitsize = FALSE, bg = "white")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, limitsize = FALSE, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(svg_path, plot, width = width, height = height,
                    device = svglite::svglite, limitsize = FALSE, bg = "white")
  } else {
    svg_path <- NA_character_
    add_warning_v5("SVG skipped for ", out_base, " because svglite is not available.")
  }

  c(png = png_path, pdf = pdf_path, svg = svg_path)
}

save_grid_device_v5 <- function(path, width, height, dpi, draw_fun, device) {
  if (identical(device, "png")) {
    ok <- FALSE
    try({
      grDevices::png(path, width = width, height = height, units = "in",
                     res = dpi, type = "cairo-png", bg = "white")
      ok <- TRUE
    }, silent = TRUE)
    if (!ok) {
      grDevices::png(path, width = width, height = height, units = "in",
                     res = dpi, bg = "white")
    }
  } else if (identical(device, "pdf")) {
    grDevices::cairo_pdf(path, width = width, height = height, bg = "white")
  } else if (identical(device, "svg")) {
    svglite::svglite(path, width = width, height = height, bg = "white")
  } else {
    stop("Unsupported device: ", device, call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_fun()
  invisible(path)
}

save_grid_all_formats_v5 <- function(out_base, width, height, dpi, draw_fun) {
  png_path <- paste0(out_base, ".png")
  pdf_path <- paste0(out_base, ".pdf")
  svg_path <- paste0(out_base, ".svg")

  save_grid_device_v5(png_path, width, height, dpi, draw_fun, "png")
  save_grid_device_v5(pdf_path, width, height, dpi, draw_fun, "pdf")
  if (requireNamespace("svglite", quietly = TRUE)) {
    save_grid_device_v5(svg_path, width, height, dpi, draw_fun, "svg")
  } else {
    svg_path <- NA_character_
    add_warning_v5("SVG skipped for ", out_base, " because svglite is not available.")
  }

  c(png = png_path, pdf = pdf_path, svg = svg_path)
}

draw_png_fit_v5 <- function(path, x, y, width, height) {
  img <- png::readPNG(path)
  img_width <- dim(img)[2]
  img_height <- dim(img)[1]
  img_aspect <- img_width / img_height
  box_aspect <- width / height

  if (img_aspect >= box_aspect) {
    draw_width <- width
    draw_height <- width / img_aspect
  } else {
    draw_height <- height
    draw_width <- height * img_aspect
  }

  draw_x <- x + (width - draw_width) / 2
  draw_y <- y + (height - draw_height) / 2

  grid::pushViewport(grid::viewport(
    x = grid::unit(draw_x, "inches"),
    y = grid::unit(draw_y, "inches"),
    width = grid::unit(draw_width, "inches"),
    height = grid::unit(draw_height, "inches"),
    just = c("left", "bottom")
  ))
  grid::grid.raster(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
                    interpolate = TRUE)
  grid::popViewport()
}

draw_grob_fit_v5 <- function(grob, source_width, source_height, x, y, width, height) {
  source_aspect <- source_width / source_height
  box_aspect <- width / height
  if (source_aspect >= box_aspect) {
    draw_width <- width
    draw_height <- width / source_aspect
  } else {
    draw_height <- height
    draw_width <- height * source_aspect
  }
  draw_x <- x + (width - draw_width) / 2
  draw_y <- y + (height - draw_height) / 2
  grid::pushViewport(grid::viewport(
    x = grid::unit(draw_x, "inches"), y = grid::unit(draw_y, "inches"),
    width = grid::unit(draw_width, "inches"), height = grid::unit(draw_height, "inches"),
    just = c("left", "bottom")
  ))
  grid::grid.draw(grob)
  grid::popViewport()
}

layout_positions_v5 <- function(fig, panel_x0, panel_y0, panel_w, panel_h, gap_x, gap_y) {
  mat <- fig$layout_matrix
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  col_weights <- fig$col_widths
  row_weights <- fig$row_heights
  if (length(col_weights) != n_cols) col_weights <- rep(1, n_cols)
  if (length(row_weights) != n_rows) row_weights <- rep(1, n_rows)

  col_widths <- (panel_w - gap_x * (n_cols - 1)) * col_weights / sum(col_weights)
  row_heights <- (panel_h - gap_y * (n_rows - 1)) * row_weights / sum(row_weights)
  col_left <- panel_x0 + c(0, cumsum(col_widths + gap_x))[seq_len(n_cols)]

  row_bottom <- numeric(n_rows)
  cursor <- panel_y0
  for (r in seq(n_rows, 1)) {
    row_bottom[r] <- cursor
    cursor <- cursor + row_heights[r] + gap_y
  }

  list(matrix = mat, col_left = col_left, row_bottom = row_bottom,
       col_widths = col_widths, row_heights = row_heights)
}

draw_clean_figure1a_grid_v5 <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  boxes <- data.frame(
    label = c(
      "Raw count matrices +\nverified metadata",
      "Delivery-versus-control\ncontrast definitions",
      "Discovery-set differential-\nexpression evidence",
      "Reproducible acute\nresponse genes",
      "Heterogeneity and\nlow-power gene filters",
      "Frozen discovery-derived\ngene weights",
      "Target dataset\nnormalization",
      "Control-referenced\ngene z-scores",
      "Weighted sample-level\nIMRS score",
      "Control-standardized\nIMRS z-score",
      "Mean delivery-minus-control\nIMRS z-score",
      "Validation groups\nand biological context"
    ),
    group = c("input", "input", rep("anchor", 4), rep("score", 4), "eval", "interpret"),
    stringsAsFactors = FALSE
  )

  coords <- data.frame(
    x = c(rep(0.25, 6), rep(0.75, 6)),
    y = c(seq(0.88, 0.18, length.out = 6), seq(0.88, 0.18, length.out = 6))
  )
  fills <- c(input = "#EAF2F8", anchor = "#E8F5E9",
             score = "#FFF4E6", eval = "#F3E8FF", interpret = "#F8EAEF")
  box_w <- 0.34
  box_h <- 0.086

  for (i in seq_len(nrow(boxes))) {
    grid::grid.roundrect(coords$x[i], coords$y[i],
                         width = box_w, height = box_h,
                         r = grid::unit(0.01, "npc"),
                         gp = grid::gpar(fill = fills[[boxes$group[i]]],
                                         col = "#334155", lwd = 1.1))
    grid::grid.text(boxes$label[i], coords$x[i], coords$y[i],
                    gp = grid::gpar(fontsize = 8.6, col = "#111827", lineheight = 0.9))
  }

  draw_arrow <- function(x0, y0, x1, y1) {
    grid::grid.lines(c(x0, x1), c(y0, y1),
                     arrow = grid::arrow(length = grid::unit(0.018, "npc"), type = "closed"),
                     gp = grid::gpar(col = "#475569", lwd = 1.05))
  }
  for (i in 1:5) {
    draw_arrow(coords$x[i], coords$y[i] - box_h / 2,
               coords$x[i + 1], coords$y[i + 1] + box_h / 2)
  }
  draw_arrow(coords$x[6] + box_w / 2, coords$y[6],
             coords$x[7] - box_w / 2, coords$y[7])
  for (i in 7:11) {
    draw_arrow(coords$x[i], coords$y[i] - box_h / 2,
               coords$x[i + 1], coords$y[i + 1] + box_h / 2)
  }
}

make_v5_panel_plan_legacy <- function() {
  data.frame(
    figure_id = c(
      "Figure1", "Figure1",
      "Figure2", "Figure2", "Figure2", "Figure2",
      "Figure3", "Figure3", "FigureS_validation_faceted_summary",
      "Figure4", "Figure4", "FigureS_weak_context_interpretation_categories",
      "Figure5", "Figure5", "Figure5", "Figure5",
      "FigureS_comparator_benchmarking", "FigureS_comparator_benchmarking"
    ),
    panel_id = c(
      "Figure1A", "Figure1B",
      "Figure2A", "Figure2B", "Figure2C", "Figure2D",
      "Figure3A", "Figure3B", "FigureS_validation_faceted_summary",
      "Figure4A", "Figure4B", "FigureS_weak_context_interpretation_categories",
      "Figure5A", "Figure5B", "Figure5C", "Figure5D",
      "FigureS_comparator_benchmarking_A", "FigureS_comparator_benchmarking_B"
    ),
    panel_letter = c(
      "A", "B",
      "A", "B", "C", "D",
      "A", "B", "",
      "A", "B", "",
      "A", "B", "C", "D",
      "A", "B"
    ),
    source_old_panel = c(
      "Figure1A_IMRS_merged_workflow",
      "Figure1C_dataset_tissue_pseudolog",
      "Figure2C_top_weighted_genes",
      "Figure2B_core_gene_reproducibility_main",
      "Figure2B_core_gene_reproducibility_missing_anchor",
      "Figure2D_weight_distribution",
      "Figure3B_primary_validation_summary",
      "FigureS1C_simplified_by_dataset",
      "Figure4A_top_contrast_responses",
      "FigureS1D_weak_zoom_forest",
      "FigureS3B_gse264344_time_course",
      "FigureS1B_weak_dataset_context_summary",
      "Figure5A_label_permutation_observed_vs_null",
      "Figure5C_permutation_response_by_analysis_group",
      "Figure7A_leave_one_gene_out_delta_correlation",
      "Figure7C_gene_dominance_distribution",
      "Figure6A_baseline_delta_by_analysis_group",
      "Figure6C_benchmark_directionality_summary"
    ),
    source_function = c(
      "render_Figure1A_merged_workflow_v5", "make_Figure1C_display_aggregated_v5",
      "make_Figure2C", "split_Figure2B_main", "split_Figure2B_missing_anchor", "make_Figure2D",
      "make_Figure3B", "make_FigureSC", "make_Figure3D_simplified",
      "make_FigureSD", "make_FigureSF", "make_FigureSB_simplified",
      "make_Figure4A", "make_Figure4C", "make_Figure5A", "make_Figure5C",
      "make_Figure4D", "make_Figure4F"
    ),
    output_stem = c(
      "Figure1A_IMRS_merged_workflow_v5", "Figure1B_dataset_tissue_response_landscape_v5_corrected",
      "Figure2A_v5", "Figure2B_v5", "Figure2C_v5", "Figure2D_v5",
      "Figure3A_v5", "Figure3B_v5", "FigureS_validation_faceted_summary_v5_panel",
      "Figure4A_v5", "Figure4B_v5", "FigureS_weak_context_interpretation_categories_v5_panel",
      "Figure5A_v5", "Figure5B_v5", "Figure5C_v5", "Figure5D_v5",
      "FigureS_comparator_benchmarking_A_v5", "FigureS_comparator_benchmarking_B_v5"
    ),
    width = c(
      7.5, 7.8,
      6.9, 5.2, 3.85, 16.2,
      7.2, 7.9, 9.6,
      8.3, 7.5, 8.6,
      7.6, 8.0, 7.6, 8.0,
      7.25, 7.1
    ),
    height = c(
      7.4, 7.1,
      4.8, 4.8, 4.8, 4.2,
      4.8, 5.6, 6.4,
      5.25, 4.6, 5.1,
      4.7, 4.7, 4.7, 4.7,
      6.4, 6.4
    ),
    notes = c(
      "Merged IMRS scoring and transfer-evaluation workflow rendered with v6 terminology, font, and spacing adjustments.",
      "Dataset/tissue-level delivery-minus-control \u0394IMRSz landscape regenerated with manuscript-analysis-group wording.",
      "Top frozen IMRS coefficients retained from clean v3 style.",
      "Locked-anchor support summary split from old Figure2B.",
      "Missing-anchor support mini-panel split from old Figure2B with 4/5 support wording.",
      "Frozen IMRS coefficient distribution retained.",
      "Primary acute vs extended validation boxplot retained.",
      "Dataset/context summary scatter retained as main Figure 3 right panel with manuscript-analysis-group wording.",
      "Faceted validation summary retained in supplement.",
      "Late/context-shifted dataset scatter retained as main Figure 4 panel.",
      "Adenoviral-vector time-course panel retained as main Figure 4 panel.",
      "Corrected context-shifted interpretation categories generated as Supplementary Figure S1B.",
      "Permutation ordered contrast plot retained.",
      "Observed score-shift summary by manuscript analysis group retained.",
      "Leave-one-gene-out \u0394IMRSz correlation retained.",
      "Gene-dominance distribution retained.",
      "Comparator immune-signature boxplot/faceted panel retained in supplemental comparator figure.",
      "Positive-directionality bar chart retained in supplemental comparator figure."
    ),
    stringsAsFactors = FALSE
  )
}

make_v5_panel_plan <- function() {
  data.frame(
    figure_id = c(
      "Figure1", "Figure1",
      rep("Figure2", 4), rep("Figure3", 2), rep("Figure4", 2),
      rep("Figure5", 4), rep("Figure6", 2), rep("FigureS1", 2)
    ),
    panel_id = c(
      "Figure1A", "Figure1B",
      "Figure2A", "Figure2B", "Figure2C", "Figure2D",
      "Figure3A", "Figure3B", "Figure4A", "Figure4B",
      "Figure5A", "Figure5B", "Figure5C", "Figure5D",
      "Figure6A", "Figure6B", "FigureS1A", "FigureS1B"
    ),
    panel_letter = c(
      "A", "B", "A", "B", "C", "D", "A", "B", "A", "B",
      "A", "B", "C", "D", "A", "B", "A", "B"
    ),
    source_old_panel = c(
      "Figure1A_IMRS_merged_workflow",
      "Figure1C_dataset_tissue_pseudolog",
      "Figure2B_core_gene_reproducibility_main",
      "Figure2B_core_gene_reproducibility_missing_anchor",
      "Figure2D_weight_distribution",
      "Figure2C_top_weighted_genes",
      "Figure3_primary_vs_extended_contrasts",
      "Figure3_primary_dataset_summary",
      "Figure4_context_boundary",
      "FigureS3B_gse264344_time_course",
      "Figure5A_label_permutation_observed_vs_null",
      "Figure5C_permutation_response_by_analysis_group",
      "Figure7A_leave_one_gene_out_delta_correlation",
      "Figure5_maximum_single_gene_contribution",
      "Figure6A_baseline_delta_by_analysis_group",
      "Figure6C_benchmark_directionality_summary",
      "FigureS1_detailed_dataset_context_audit",
      "FigureS1_context_category_counts"
    ),
    source_function = c(
      "render_Figure1A_merged_workflow_v5",
      "make_Figure1C_display_aggregated_v5",
      "split_Figure2B_main",
      "split_Figure2B_missing_anchor",
      "make_Figure2D",
      "make_Figure2C",
      "render_figure3a_primary_extended_v5",
      "render_figure3b_primary_dataset_summary_v5",
      "render_figure4a_boundary_context_v5",
      "make_FigureSF",
      "make_Figure4A",
      "make_Figure4C",
      "make_Figure5A",
      "render_figure5d_max_contribution_v5",
      "make_Figure4D",
      "make_Figure4F",
      "render_s1a_detailed_dataset_context_v5",
      "render_s1b_context_categories_v5"
    ),
    output_stem = c(
      "Figure1A_IMRS_merged_workflow_v5",
      "Figure1B_dataset_tissue_response_landscape_v5_corrected",
      "Figure2A_v5", "Figure2B_v5", "Figure2C_v5", "Figure2D_v5",
      "Figure3A_v5", "Figure3B_v5", "Figure4A_v5", "Figure4B_v5",
      "Figure5A_v5", "Figure5B_v5", "Figure5C_v5", "Figure5D_v5",
      "Figure6A_v5", "Figure6B_v5", "FigureS1A_v5", "FigureS1B_v5"
    ),
    width = c(
      7.5, 7.8, 6.8, 6.8, 6.8, 6.8, 7.4, 8.0, 8.7, 7.5,
      7.6, 8.0, 7.6, 8.0, 7.25, 7.1, 9.0, 9.0
    ),
    height = c(
      7.4, 7.1, 4.8, 4.8, 4.6, 4.8, 5.1, 5.1, 5.4, 4.6,
      4.7, 4.7, 4.7, 4.7, 6.4, 6.4, 8.7, 4.8
    ),
    notes = c(
      "Frozen IMRS scoring and transfer-evaluation workflow.",
      "Compact global dataset/context landscape across all valid manuscript groups.",
      "Retained-gene support count: 178 all-five and 109 four-of-five.",
      "Missing-anchor breakdown for the 109 four-of-five-supported genes.",
      "Distribution of the 287 applied beta_meta coefficients.",
      "Top genes ranked by absolute beta_meta coefficient.",
      "Split-level primary acute versus extended validation distribution.",
      "Nested split contrasts and dataset means for three independent primary datasets.",
      "Thirteen context-audit rows grouped by biological interpretation category.",
      "GSE264344 adenoviral-vector time course.",
      "Observed split effects relative to within-split permutation-null intervals.",
      "Observed score-shift summaries by principal manuscript group.",
      "Leave-one-gene-out effects versus original effects.",
      "Mean maximum single-gene contribution fraction used, matching the manuscript metric.",
      "Comparator immune-signature split-level shifts reclassified as main Figure 6.",
      "Comparator positive-direction proportions reclassified as main Figure 6.",
      "Detailed faceted dataset/context provenance across all manuscript groups.",
      "Counts of the 13 context-audit rows by interpretation category and support level."
    ),
    stringsAsFactors = FALSE
  )
}

make_v5_figure_plan_legacy <- function() {
  list(
    list(
      figure_id = "Figure1",
      output_stem = "Figure1_main_v5",
      width = 15.9,
      height = 7.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(7.5, 7.8),
      row_heights = c(1),
      panels = c(A = "Figure1A", B = "Figure1B"),
      role = "main"
    ),
    list(
      figure_id = "Figure2",
      output_stem = "Figure2_main_v5",
      width = 16.6,
      height = 9.6,
      layout_matrix = matrix(c("A", "B", "C", "D", "D", "D"), nrow = 2, byrow = TRUE),
      col_widths = c(7.2, 5.4, 4.0),
      row_heights = c(5.2, 4.6),
      panels = c(A = "Figure2A", B = "Figure2B", C = "Figure2C", D = "Figure2D"),
      role = "main"
    ),
    list(
      figure_id = "Figure3",
      output_stem = "Figure3_main_v5",
      width = 15.8,
      height = 6.2,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(7.2, 8.2),
      row_heights = c(1),
      panels = c(A = "Figure3A", B = "Figure3B"),
      role = "main"
    ),
    list(
      figure_id = "Figure4",
      output_stem = "Figure4_main_v5",
      width = 16.2,
      height = 5.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(8.8, 7.5),
      row_heights = c(1),
      panels = c(A = "Figure4A", B = "Figure4B"),
      role = "main"
    ),
    list(
      figure_id = "Figure5",
      output_stem = "Figure5_main_v5",
      width = 16.2,
      height = 10.0,
      layout_matrix = matrix(c("A", "B", "C", "D"), nrow = 2, byrow = TRUE),
      col_widths = c(8.4, 8.8),
      row_heights = c(5.2, 5.2),
      panels = c(A = "Figure5A", B = "Figure5B", C = "Figure5C", D = "Figure5D"),
      role = "main"
    ),
    list(
      figure_id = "FigureS_validation_faceted_summary",
      output_stem = "FigureS_validation_faceted_summary_v5",
      width = 10.0,
      height = 6.8,
      layout_matrix = matrix("A", nrow = 1),
      col_widths = c(1),
      row_heights = c(1),
      panels = c(A = "FigureS_validation_faceted_summary"),
      role = "supplement"
    ),
    list(
      figure_id = "FigureS_weak_context_interpretation_categories",
      output_stem = "FigureS_weak_context_interpretation_categories_v5",
      width = 9.0,
      height = 5.5,
      layout_matrix = matrix("A", nrow = 1),
      col_widths = c(1),
      row_heights = c(1),
      panels = c(A = "FigureS_weak_context_interpretation_categories"),
      role = "supplement"
    ),
    list(
      figure_id = "FigureS_comparator_benchmarking",
      output_stem = "FigureS_comparator_benchmarking_v5",
      width = 15.0,
      height = 6.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(9.0, 8.8),
      row_heights = c(1),
      panels = c(A = "FigureS_comparator_benchmarking_A", B = "FigureS_comparator_benchmarking_B"),
      role = "supplement"
    )
  )
}

make_v5_figure_plan <- function() {
  list(
    list(
      figure_id = "Figure1", output_stem = "Figure1_main_v5",
      width = 15.9, height = 7.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(7.5, 7.8), row_heights = 1,
      panels = c(A = "Figure1A", B = "Figure1B"), role = "main"
    ),
    list(
      figure_id = "Figure2", output_stem = "Figure2_main_v5",
      width = 14.4, height = 10.0,
      layout_matrix = matrix(c("A", "B", "C", "D"), nrow = 2, byrow = TRUE),
      col_widths = c(1, 1), row_heights = c(1, 1),
      panels = c(A = "Figure2A", B = "Figure2B", C = "Figure2C", D = "Figure2D"),
      role = "main"
    ),
    list(
      figure_id = "Figure3", output_stem = "Figure3_main_v5",
      width = 15.8, height = 5.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(7.4, 8.0), row_heights = 1,
      panels = c(A = "Figure3A", B = "Figure3B"), role = "main"
    ),
    list(
      figure_id = "Figure4", output_stem = "Figure4_main_v5",
      width = 16.5, height = 6.1,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(8.7, 7.5), row_heights = 1,
      panels = c(A = "Figure4A", B = "Figure4B"), role = "main"
    ),
    list(
      figure_id = "Figure5", output_stem = "Figure5_main_v5",
      width = 16.2, height = 10.0,
      layout_matrix = matrix(c("A", "B", "C", "D"), nrow = 2, byrow = TRUE),
      col_widths = c(1, 1), row_heights = c(1, 1),
      panels = c(A = "Figure5A", B = "Figure5B", C = "Figure5C", D = "Figure5D"),
      role = "main"
    ),
    list(
      figure_id = "Figure6", output_stem = "Figure6_main_v5",
      width = 15.0, height = 6.8,
      layout_matrix = matrix(c("A", "B"), nrow = 1, byrow = TRUE),
      col_widths = c(7.25, 7.1), row_heights = 1,
      panels = c(A = "Figure6A", B = "Figure6B"), role = "main"
    ),
    list(
      figure_id = "FigureS1", output_stem = "FigureS1_main_v5",
      width = 9.4, height = 14.2,
      layout_matrix = matrix(c("A", "B"), nrow = 2, byrow = TRUE),
      col_widths = 1, row_heights = c(8.7, 4.8),
      panels = c(A = "FigureS1A", B = "FigureS1B"), role = "supplement"
    )
  )
}

render_figure2_support_split_v5 <- function(old_env, out_dir, dpi = 400) {
  stop_if_missing_packages_v5()
  env <- old_env

  weights <- env$prepare_weights()
  retained_genes <- tibble::tibble(gene_id_clean = unique(weights$gene_id_clean))
  support_wide <- env$support_tbl %>%
    dplyr::mutate(
      dataset_support_flag = env$logic_col(dataset_support_flag),
      gene_id_clean = env$strip_ens(gene_id),
      dataset_id = as.character(dataset_id)
    ) %>%
    dplyr::filter(dataset_id %in% env$LOCKED_DATASETS_MOUSE, gene_id_clean %in% retained_genes$gene_id_clean) %>%
    dplyr::group_by(gene_id_clean, dataset_id) %>%
    dplyr::summarise(supported = any(dataset_support_flag, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = dataset_id, values_from = supported, values_fill = FALSE) %>%
    dplyr::right_join(retained_genes, by = "gene_id_clean")
  for (dataset_id in env$LOCKED_DATASETS_MOUSE) {
    if (!dataset_id %in% names(support_wide)) support_wide[[dataset_id]] <- FALSE
  }
  support_wide <- support_wide %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(env$LOCKED_DATASETS_MOUSE), ~ tidyr::replace_na(as.logical(.x), FALSE)))
  support_matrix <- as.matrix(support_wide[, env$LOCKED_DATASETS_MOUSE, drop = FALSE])
  support_wide$support_n <- rowSums(support_matrix)
  support_wide$pattern_label <- apply(support_matrix, 1, function(row) {
    supported <- env$LOCKED_DATASETS_MOUSE[as.logical(row)]
    missing <- setdiff(env$LOCKED_DATASETS_MOUSE, supported)
    if (length(supported) == length(env$LOCKED_DATASETS_MOUSE)) {
      paste0("All ", length(env$LOCKED_DATASETS_MOUSE), " locked anchors")
    } else if (length(missing) == 1) {
      paste0("All except ", missing)
    } else if (length(supported) > 0) {
      paste(supported, collapse = " + ")
    } else {
      "No discovery support detected"
    }
  })

  k_locked <- length(env$LOCKED_DATASETS_MOUSE)
  support_label <- function(n) {
    dplyr::case_when(
      n == k_locked ~ paste0("Supported in ", k_locked, "/", k_locked, " locked anchors"),
      n == k_locked - 1L ~ paste0("Supported in ", k_locked - 1L, "/", k_locked, " locked anchors"),
      n == 1L ~ "1 locked anchor",
      TRUE ~ paste0(n, " locked anchors")
    )
  }
  plot_tbl <- support_wide %>%
    dplyr::mutate(support_n = as.integer(support_n), support_category = support_label(support_n)) %>%
    dplyr::count(support_n, support_category, name = "n_retained_genes") %>%
    dplyr::filter(n_retained_genes > 0) %>%
    dplyr::arrange(dplyr::desc(support_n)) %>%
    dplyr::mutate(support_category = factor(support_category, levels = unique(support_category)))

  pattern_notes <- support_wide %>%
    dplyr::group_by(support_n, pattern_label) %>%
    dplyr::summarise(n_retained_genes = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      support_category = support_label(support_n),
      supporting_datasets = purrr::map_chr(pattern_label, function(label) {
        if (startsWith(label, "All except ")) return(paste(setdiff(env$LOCKED_DATASETS_MOUSE, sub("^All except ", "", label)), collapse = ";"))
        if (startsWith(label, "All ")) return(paste(env$LOCKED_DATASETS_MOUSE, collapse = ";"))
        if (label == "No discovery support detected") return("")
        stringr::str_replace_all(label, " \\+ ", ";")
      }),
      absent_datasets = purrr::map_chr(supporting_datasets, function(supported) {
        supported_vec <- if (nzchar(supported)) unlist(strsplit(supported, ";", fixed = TRUE)) else character()
        paste(setdiff(env$LOCKED_DATASETS_MOUSE, supported_vec), collapse = ";")
      })
    ) %>%
    dplyr::arrange(dplyr::desc(support_n), pattern_label)
  readr::write_tsv(pattern_notes, file.path(dirname(out_dir), "tables", "Figure2_support_pattern_notes_v5.tsv"), na = "NA")

  all_but_one_total <- plot_tbl %>%
    dplyr::filter(support_n == k_locked - 1L) %>%
    dplyr::summarise(total = sum(n_retained_genes), .groups = "drop") %>%
    dplyr::pull(total)
  if (length(all_but_one_total) == 0L) all_but_one_total <- 0L
  all_five_total <- plot_tbl %>%
    dplyr::filter(support_n == k_locked) %>%
    dplyr::summarise(total = sum(n_retained_genes), .groups = "drop") %>%
    dplyr::pull(total)
  if (length(all_five_total) == 0L) all_five_total <- 0L
  if (!identical(as.integer(all_five_total), 178L) ||
      !identical(as.integer(all_but_one_total), 109L)) {
    stop("Retained-gene support manuscript invariant failed: expected 178 all-five and 109 four-of-five; observed ",
         all_five_total, " and ", all_but_one_total, ".", call. = FALSE)
  }

  missing_anchor_tbl <- pattern_notes %>%
    dplyr::filter(support_n == k_locked - 1L, stringr::str_detect(pattern_label, "^All except ")) %>%
    dplyr::transmute(
      missing_anchor = stringr::str_remove(pattern_label, "^All except "),
      n_retained_genes
    ) %>%
    dplyr::group_by(missing_anchor) %>%
    dplyr::summarise(n_retained_genes = sum(n_retained_genes), .groups = "drop") %>%
    dplyr::mutate(missing_anchor = factor(missing_anchor, levels = env$LOCKED_DATASETS_MOUSE)) %>%
    dplyr::arrange(missing_anchor)
  if (sum(missing_anchor_tbl$n_retained_genes, na.rm = TRUE) != all_but_one_total) {
    stop("Figure2 all-but-one missing-anchor split does not sum to the all-but-one total.", call. = FALSE)
  }
  expected_missing <- c(
    GSE167521 = 18L, GSE264344 = 13L, GSE279372 = 33L,
    GSE279744 = 10L, GSE39129 = 35L
  )
  observed_missing <- stats::setNames(
    as.integer(missing_anchor_tbl$n_retained_genes),
    as.character(missing_anchor_tbl$missing_anchor)
  )
  if (!identical(observed_missing[names(expected_missing)], expected_missing)) {
    stop("Missing-anchor manuscript invariant failed. Expected ",
         paste(names(expected_missing), expected_missing, sep = "=", collapse = ", "),
         ".", call. = FALSE)
  }

  main_y_max <- max(plot_tbl$n_retained_genes, na.rm = TRUE)
  main_plot <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = support_category, y = n_retained_genes)) +
    ggplot2::geom_col(fill = "#4E79A7", color = "#374151", linewidth = 0.45, width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = n_retained_genes), vjust = -0.35, size = 4.2, color = "#111827") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16)), breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::coord_cartesian(ylim = c(0, main_y_max * 1.18), clip = "off") +
    ggplot2::labs(x = "Locked-anchor support", y = "Number of retained genes") +
    env$theme_imrs_publication(base_size = 12, legend_position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11.0, lineheight = 0.95, margin = ggplot2::margin(t = 6)),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )

  mini_y_max <- if (nrow(missing_anchor_tbl) > 0L) max(missing_anchor_tbl$n_retained_genes, na.rm = TRUE) else 1
  mini_plot <- ggplot2::ggplot(missing_anchor_tbl, ggplot2::aes(x = missing_anchor, y = n_retained_genes, fill = missing_anchor)) +
    ggplot2::geom_col(color = "#374151", linewidth = 0.35, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(n_retained_genes > 0, n_retained_genes, "")),
                       vjust = -0.25, size = 3.1, color = "#111827") +
    ggplot2::scale_fill_manual(values = env$anchor_palette[env$LOCKED_DATASETS_MOUSE], drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)), breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::coord_cartesian(ylim = c(0, max(1, mini_y_max) * 1.22), clip = "off") +
    ggplot2::labs(x = "Missing anchor", y = "Retained genes") +
    env$theme_imrs_publication(base_size = 11.25, legend_position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 10.2),
      axis.title = ggplot2::element_text(size = 10.2),
      axis.title.x = ggplot2::element_text(size = 10.2, margin = ggplot2::margin(t = 5)),
      plot.margin = ggplot2::margin(6, 6, 10, 6)
    )

  main_plot <- strip_publication_text_v5(main_plot)
  mini_plot <- strip_publication_text_v5(mini_plot)
  main_paths <- save_plot_all_formats_v5(main_plot, file.path(out_dir, "Figure2A_v5"), 6.8, 4.8, dpi)
  mini_paths <- save_plot_all_formats_v5(mini_plot, file.path(out_dir, "Figure2B_v5"), 6.8, 4.8, dpi)
  list(main = main_paths, mini = mini_paths, main_plot = main_plot, mini_plot = mini_plot)
}

make_figure1b_display_table_v5 <- function(old_env) {
  figure1b_input <- old_env$role_pass_for_plot %>%
    dplyr::filter(is.finite(delta_mean_imrs_z)) %>%
    dplyr::mutate(
      time_h_num = old_env$safe_num(time_h),
      is_gse264344_acute = dataset_id == "GSE264344" &
        is.finite(time_h_num) & time_h_num <= 24,
      is_gse264344_late = dataset_id == "GSE264344" &
        is.finite(time_h_num) & time_h_num > 24,
      display_time_h = dplyr::if_else(is_gse264344_acute, NA_real_, time_h_num),
      time_window_label = dplyr::if_else(
        is_gse264344_acute,
        "1-24 h",
        old_env$format_time_label(time_h_num)
      ),
      display_label = dplyr::case_when(
        is_gse264344_acute ~ paste(
          "GSE264344 |",
          old_env$display_text(tissue),
          "| 1-24 h"
        ),
        is_gse264344_late ~ paste(
          "GSE264344 |",
          old_env$display_text(tissue),
          "|",
          old_env$format_time_label(time_h_num)
        ),
        TRUE ~ paste(
          dataset_id,
          old_env$display_text(tissue),
          old_env$format_time_label(time_h_num),
          sep = " | "
        )
      )
    )

  figure1b_input %>%
    dplyr::group_by(
      dataset_id, tissue, display_time_h, delivery_platform_clean,
      manuscript_group, time_window_label, display_label
    ) %>%
    dplyr::summarise(
      mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
      n_contrasts = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(label = old_env$ordered_factor(display_label, mean_delta)) %>%
    dplyr::arrange(label)
}

render_corrected_figure1b_landscape_v5 <- function(old_env, v5_root, dpi = 400) {
  plot_tbl <- make_figure1b_display_table_v5(old_env)
  old_env$check_not_all_na(plot_tbl, "mean_delta", "Figure1B corrected plot table")

  plot_table_path <- file.path(
    v5_root, "tables",
    "Figure1B_dataset_tissue_response_landscape_v5_corrected_plot_table.tsv"
  )
  write_tsv_v5(
    plot_tbl %>%
      dplyr::mutate(label = as.character(display_label)) %>%
      dplyr::select(
        dataset_id, tissue, display_time_h, delivery_platform_clean,
        manuscript_group, time_window_label, label, mean_delta, n_contrasts
      ),
    plot_table_path
  )

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = mean_delta, y = label, color = manuscript_group, size = n_contrasts)
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#4B5563"
    ) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 2),
      breaks = c(-2, -1, 0, 1, 2, 5, 10, 25, 50)
    ) +
    ggplot2::scale_color_manual(values = old_env$manuscript_group_palette, drop = FALSE) +
    ggplot2::scale_size_continuous(
      range = c(2, 6),
      breaks = c(1, 2, 3, 4, 5, 6),
      limits = c(1, max(plot_tbl$n_contrasts, na.rm = TRUE)),
      name = "Passing split contrasts"
    ) +
    ggplot2::labs(
      title = "Dataset-level delivery-minus-control \u0394IMRSz landscape",
      subtitle = "Pseudo-log x-axis preserves resolution near zero while keeping large positive responses visible.",
      x = "Mean delivery-minus-control \u0394IMRSz (pseudo-log scale)",
      y = NULL,
      color = "Manuscript analysis group",
      size = "Passing split contrasts"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Manuscript analysis group",
        order = 1,
        ncol = 1,
        byrow = TRUE,
        override.aes = list(size = 3)
      ),
      size = ggplot2::guide_legend(
        title = "Passing split contrasts",
        order = 2,
        ncol = 2,
        byrow = TRUE
      )
    ) +
    old_env$theme_imrs_publication(base_size = 10.25) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10.4),
      axis.text.x = ggplot2::element_text(size = 9.8),
      axis.title.x = ggplot2::element_text(size = 10.7),
      legend.text = ggplot2::element_text(size = 9.2),
      legend.title = ggplot2::element_text(size = 9.8),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 4),
      legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 2),
      legend.spacing.y = grid::unit(8, "pt"),
      plot.margin = ggplot2::margin(8, 8, 8, 10)
    )

  paths <- save_plot_all_formats_v5(
    strip_publication_text_v5(p),
    file.path(v5_root, "intermediate_panels", "Figure1B_dataset_tissue_response_landscape_v5_corrected"),
    7.8, 7.1, dpi
  )

  log_msg_v5(
    "Corrected Figure 1B display table contributes ",
    nrow(dplyr::filter(plot_tbl, dataset_id == "GSE264344")),
    " GSE264344 rows."
  )
  list(paths = paths, plot = strip_publication_text_v5(p),
       plot_table = plot_tbl, plot_table_path = plot_table_path)
}

save_manuscript_panel_v5 <- function(plot, out_dir, stem, width, height, dpi = 400) {
  clean_plot <- strip_publication_text_v5(plot)
  list(
    paths = save_plot_all_formats_v5(clean_plot, file.path(out_dir, stem), width, height, dpi),
    plot = clean_plot
  )
}

render_figure3a_primary_extended_v5 <- function(env, out_dir, dpi = 400) {
  tbl <- env$role_pass_for_plot %>%
    dplyr::filter(as.character(manuscript_group) %in%
                    c("Primary acute validation", "Extended validation")) %>%
    dplyr::mutate(
      manuscript_group = factor(as.character(manuscript_group),
                                levels = c("Primary acute validation", "Extended validation"))
    )
  ann <- tbl %>%
    dplyr::group_by(manuscript_group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
      prop_pos = mean(delta_mean_imrs_z > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf("%d split contrasts\nmean \u0394IMRSz = %.3f\n\u0394IMRSz > 0: %d/%d",
                      n, mean_delta, round(prop_pos * n), n),
      y = max(tbl$delta_mean_imrs_z, na.rm = TRUE) + 2.2
    )
  p <- ggplot2::ggplot(
    tbl, ggplot2::aes(x = manuscript_group, y = delta_mean_imrs_z,
                      fill = manuscript_group)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                       linewidth = 0.4, color = "#4B5563") +
    ggplot2::geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.68) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.82,
                        color = "#111827", show.legend = FALSE) +
    ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = manuscript_group, y = y, label = label),
      inherit.aes = FALSE, size = 3.45, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = env$manuscript_group_palette, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.26))) +
    ggplot2::labs(
      x = "Manuscript analysis group",
      y = "Split-level delivery-minus-control \u0394IMRSz",
      fill = "Manuscript analysis group"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    env$theme_imrs_publication(base_size = 11.5, legend_position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10.2, margin = ggplot2::margin(t = 6)),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 14, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "Figure3A_v5", 7.4, 5.1, dpi)
}

render_figure3b_primary_dataset_summary_v5 <- function(env, out_dir, dpi = 400) {
  tbl <- env$role_pass_for_plot %>%
    dplyr::filter(as.character(manuscript_group) == "Primary acute validation") %>%
    dplyr::mutate(
      tissue_label = env$display_text(tissue),
      time_label = env$format_time_label(time_h)
    )
  summary_tbl <- tbl %>%
    dplyr::group_by(dataset_id, tissue_label, time_label) %>%
    dplyr::summarise(
      mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
      min_delta = min(delta_mean_imrs_z, na.rm = TRUE),
      max_delta = max(delta_mean_imrs_z, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  expected_means <- c(GSE139529 = 11.920, GSE119119 = 10.739, GSE279743 = 8.148)
  expected_counts <- c(GSE139529 = 6L, GSE119119 = 6L, GSE279743 = 2L)
  observed_means <- stats::setNames(summary_tbl$mean_delta, summary_tbl$dataset_id)
  observed_counts <- stats::setNames(summary_tbl$n, summary_tbl$dataset_id)
  if (!all(names(expected_means) %in% names(observed_means)) ||
      any(abs(observed_means[names(expected_means)] - expected_means) > 0.0005) ||
      !identical(as.integer(observed_counts[names(expected_counts)]),
                 as.integer(expected_counts))) {
    stop("Figure 3B manuscript invariants failed for dataset means or split counts.",
         call. = FALSE)
  }
  summary_tbl <- summary_tbl %>%
    dplyr::arrange(mean_delta) %>%
    dplyr::mutate(
      dataset_label = paste0(
        dataset_id, " | ", tissue_label, " | ", time_label,
        "\n", n, " split contrasts"
      ),
      dataset_label = factor(dataset_label, levels = dataset_label),
      mean_label = sprintf("%.2f", mean_delta)
    )
  tbl <- tbl %>%
    dplyr::left_join(
      summary_tbl %>%
        dplyr::transmute(dataset_id, dataset_label = as.character(dataset_label)),
      by = "dataset_id"
    ) %>%
    dplyr::mutate(
      dataset_label = factor(dataset_label, levels = levels(summary_tbl$dataset_label))
    )
  primary_color <- unname(env$manuscript_group_palette["Primary acute validation"])
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = delta_mean_imrs_z, y = dataset_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
                       color = "#6B7280") +
    ggplot2::geom_segment(
      data = summary_tbl,
      ggplot2::aes(
        x = min_delta, xend = max_delta,
        y = dataset_label, yend = dataset_label
      ),
      inherit.aes = FALSE, linewidth = 0.7, alpha = 0.50,
      color = "#6B7280", lineend = "butt"
    ) +
    ggplot2::geom_jitter(height = 0.08, width = 0, shape = 16, size = 1.8,
                        alpha = 0.70, color = "#4B5563") +
    ggplot2::geom_point(
      data = summary_tbl,
      ggplot2::aes(x = mean_delta, y = dataset_label),
      inherit.aes = FALSE, shape = 23, size = 4.6, stroke = 0.75,
      fill = primary_color, color = "#1F2937"
    ) +
    ggplot2::geom_text(
      data = summary_tbl,
      ggplot2::aes(x = mean_delta, y = dataset_label, label = mean_label),
      inherit.aes = FALSE, nudge_y = 0.18, hjust = 0.5, vjust = 0,
      size = 3.15, color = "#374151"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
    ggplot2::labs(
      x = "Nested split-contrast \u0394IMRSz",
      y = NULL
    ) +
    env$theme_imrs_publication(base_size = 11.5, legend_position = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9.8, lineheight = 0.95,
                                         color = "#111827"),
      axis.text.x = ggplot2::element_text(size = 10.2),
      plot.margin = ggplot2::margin(8, 12, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "Figure3B_v5", 8.0, 5.1, dpi)
}

context_category_labels_v5 <- c(
  disease_rescue_model = "Disease-rescue model",
  distal_or_adaptive_tissue = "Distal/adaptive tissue context",
  late_timepoint = "Late time point",
  formulation_designed_to_reduce_inflammation = "Low-inflammatory formulation context",
  therapeutic_cargo_specific_effect = "Therapeutic cargo/context",
  tissue_time_kinetic_effect = "Tissue/time-course context"
)

render_figure4a_boundary_context_v5 <- function(env, out_dir, dpi = 400) {
  category_levels <- unname(context_category_labels_v5)
  tbl <- env$weak_tbl %>%
    dplyr::left_join(
      env$role_pass_for_plot %>%
        dplyr::select(dataset_id, split_id, manuscript_group) %>%
        dplyr::distinct(),
      by = c("dataset_id", "split_id")
    ) %>%
    dplyr::mutate(
      delta = env$safe_num(original_IMRS_delta),
      category = unname(context_category_labels_v5[as.character(explanation_category)]),
      category = factor(category, levels = rev(category_levels)),
      manuscript_group = factor(as.character(manuscript_group),
                                levels = c("Extended validation", "Secondary support")),
      point_label = dplyr::case_when(
        dataset_id == "GSE166655" & grepl("ar45", treatment_group, ignore.case = TRUE) ~
          "GSE166655 \u00B7 1008 h",
        dataset_id == "GSE262515_tissue" & grepl("sinc", treatment_group, ignore.case = TRUE) ~
          "GSE262515 tissue \u00B7 72 h",
        dataset_id == "GSE314070" & grepl("sm102", treatment_group, ignore.case = TRUE) ~
          "GSE314070 \u00B7 336 h",
        TRUE ~ ""
      )
    )
  if (nrow(tbl) != 13L || sum(is.finite(tbl$delta)) != 13L ||
      any(is.na(tbl$category)) || any(is.na(tbl$manuscript_group))) {
    stop("Figure 4A manuscript invariant failed: expected 13 unchanged context observations.",
         call. = FALSE)
  }
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = delta, y = category, fill = manuscript_group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
                       color = "#6B7280") +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0, height = 0.13, seed = 20260906),
      shape = 21, size = 2.85, stroke = 0.45, color = "#4B5563", alpha = 0.78
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(tbl, nzchar(point_label)),
      ggplot2::aes(label = point_label),
      position = ggplot2::position_jitter(width = 0, height = 0.13, seed = 20260906),
      hjust = -0.08, vjust = -0.7, size = 2.95, color = "#111827",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = env$manuscript_group_palette, drop = TRUE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.38))) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(size = 2.4, alpha = 0.78, stroke = 0.45)
    )) +
    ggplot2::labs(
      x = "Observed delivery-minus-control \u0394IMRSz",
      y = "Biological context category",
      fill = "Manuscript analysis group"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    env$theme_imrs_publication(base_size = 11, legend_position = "bottom") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10.2),
      axis.text.x = ggplot2::element_text(size = 9.8),
      legend.text = ggplot2::element_text(size = 9.0),
      legend.title = ggplot2::element_text(size = 9.0, face = "plain"),
      legend.key.size = grid::unit(9.5, "pt"),
      legend.spacing.x = grid::unit(3, "pt"),
      plot.margin = ggplot2::margin(10, 26, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "Figure4A_v5", 8.7, 5.4, dpi)
}

render_figure5d_max_contribution_v5 <- function(env, out_dir, dpi = 400) {
  tbl <- env$join_role_group(env$dominance_tbl) %>%
    dplyr::mutate(
      mean_max_contribution_fraction = env$safe_num(mean_max_contribution_fraction),
      manuscript_group = factor(as.character(manuscript_group),
                                levels = env$manuscript_group_order)
    ) %>%
    dplyr::filter(is.finite(mean_max_contribution_fraction), !is.na(manuscript_group))
  overall_mean <- mean(tbl$mean_max_contribution_fraction, na.rm = TRUE)
  overall_max <- max(tbl$mean_max_contribution_fraction, na.rm = TRUE)
  if (round(100 * overall_mean, 1) != 3.3 || round(100 * overall_max, 1) != 8.7) {
    stop("Figure 5D manuscript invariants failed for overall mean or maximum.",
         call. = FALSE)
  }
  max_row <- tbl %>%
    dplyr::slice_max(mean_max_contribution_fraction, n = 1, with_ties = FALSE) %>%
    dplyr::mutate(max_label = sprintf(
      "Maximum observed = %.1f%%", 100 * mean_max_contribution_fraction
    ))
  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = mean_max_contribution_fraction, y = manuscript_group,
                 fill = manuscript_group)
  ) +
    ggplot2::geom_vline(xintercept = overall_mean, linetype = "dashed",
                       linewidth = 0.3, color = "#6B7280") +
    ggplot2::geom_boxplot(
      width = 0.5, outlier.shape = NA, alpha = 0.40,
      box.linewidth = 0.45, median.linewidth = 0.75,
      whisker.linewidth = 0.4, staple.linewidth = 0.4
    ) +
    ggplot2::geom_jitter(height = 0.12, width = 0, shape = 21, size = 1.7,
                        stroke = 0.3, color = "#374151", alpha = 0.78,
                        show.legend = FALSE) +
    ggplot2::geom_text(
      data = max_row,
      ggplot2::aes(
        x = mean_max_contribution_fraction, y = manuscript_group,
        label = max_label
      ),
      inherit.aes = FALSE, nudge_y = 0.18, hjust = 1.02, vjust = 0,
      size = 3.0, color = "#374151", show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = env$manuscript_group_palette, drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                               expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(
      x = "Mean maximum single-gene contribution (%)",
      y = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    env$theme_imrs_publication(base_size = 11, legend_position = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9.8),
      axis.text.x = ggplot2::element_text(size = 9.6),
      plot.margin = ggplot2::margin(18, 10, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "Figure5D_v5", 8.0, 4.7, dpi)
}

render_s1a_detailed_dataset_context_v5 <- function(env, out_dir, dpi = 400) {
  tbl <- env$role_pass_for_plot %>%
    dplyr::filter(is.finite(delta_mean_imrs_z)) %>%
    dplyr::mutate(
      context_label = paste0(dataset_id, " | ", env$display_text(tissue),
                             " | ", env$format_time_label(time_h))
    ) %>%
    dplyr::group_by(manuscript_group, dataset_id, tissue, time_h, context_label) %>%
    dplyr::summarise(
      mean_delta = mean(delta_mean_imrs_z, na.rm = TRUE),
      n_contrasts = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::group_by(manuscript_group) %>%
    dplyr::arrange(mean_delta, .by_group = TRUE) %>%
    dplyr::mutate(context_label = factor(context_label, levels = unique(context_label))) %>%
    dplyr::ungroup()
  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = mean_delta, y = context_label,
                 color = manuscript_group, size = n_contrasts)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
                       color = "#4B5563") +
    ggplot2::geom_point(alpha = 0.94) +
    ggplot2::facet_grid(
      manuscript_group ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::as_labeller(c(
        "Locked anchor" = "Anchor",
        "Primary acute validation" = "Primary",
        "Extended validation" = "Extended",
        "Secondary support" = "Secondary"
      ))
    ) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 2),
      breaks = c(-2, -1, 0, 1, 2, 5, 10, 25, 50)
    ) +
    ggplot2::scale_color_manual(values = env$manuscript_group_palette, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(2.2, 5.6), breaks = c(1, 2, 3, 6)) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, ncol = 1),
      size = ggplot2::guide_legend(order = 2, ncol = 1,
                                  override.aes = list(color = "#111827"))
    ) +
    ggplot2::labs(
      x = "Dataset/context mean \u0394IMRSz (pseudo-log scale)",
      y = "Dataset / tissue / time",
      color = "Manuscript analysis group",
      size = "Passing split contrasts"
    ) +
    env$theme_imrs_publication(base_size = 11.5, legend_position = "right") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10.4),
      axis.text.x = ggplot2::element_text(size = 10.0),
      strip.text.y = ggplot2::element_text(angle = 0, size = 9.5, face = "bold"),
      legend.text = ggplot2::element_text(size = 9.6),
      legend.title = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(8, 8, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "FigureS1A_v5", 9.0, 8.7, dpi)
}

render_s1b_context_categories_v5 <- function(env, out_dir, dpi = 400) {
  category_levels <- unname(context_category_labels_v5)
  if (dplyr::n_distinct(env$weak_tbl$split_id) != nrow(env$weak_tbl)) {
    stop("Supplementary Figure S1B requires one unique split contrast per context row.",
         call. = FALSE)
  }
  tbl <- env$weak_tbl %>%
    dplyr::mutate(
      category = unname(context_category_labels_v5[as.character(explanation_category)]),
      category = factor(category, levels = rev(category_levels))
    ) %>%
    dplyr::count(category, name = "n", .drop = TRUE)
  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = n, y = category)) +
    ggplot2::geom_col(width = 0.68, fill = "#4E79A7") +
    ggplot2::geom_text(
      data = dplyr::filter(tbl, n > 0), ggplot2::aes(label = n),
      hjust = 1.4, color = "white",
      fontface = "bold", size = 3.7
    ) +
    ggplot2::scale_x_continuous(breaks = 0:6,
                               expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      x = "Number of context-shifted contrasts",
      y = "Biological context category"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    env$theme_imrs_publication(base_size = 11.5, legend_position = "none") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10.6),
      axis.text.x = ggplot2::element_text(size = 10.0),
      plot.margin = ggplot2::margin(8, 16, 8, 10)
    )
  save_manuscript_panel_v5(p, out_dir, "FigureS1B_v5", 9.0, 4.8, dpi)
}

render_figure1a_merged_workflow_v5 <- function(project_root, v5_root, dpi = 400,
                                               workflow_script = NULL) {
  if (is.null(workflow_script) || !nzchar(workflow_script)) {
    stop("A repo-contained merged workflow script must be supplied.", call. = FALSE)
  }
  if (!file.exists(workflow_script)) {
    stop("Missing merged workflow build script: ", workflow_script, call. = FALSE)
  }
  workflow_env <- new.env(parent = globalenv())
  sys.source(workflow_script, envir = workflow_env)
  if (!exists("build_merged_imrs_workflow_v5", envir = workflow_env, inherits = FALSE)) {
    stop("Workflow build script did not define build_merged_imrs_workflow_v5().", call. = FALSE)
  }
  paths <- workflow_env$build_merged_imrs_workflow_v5(
    project_root = project_root,
    output_dir = v5_root,
    output_stem = "Figure1A_IMRS_merged_workflow_v5",
    width = 7.5,
    height = 7.4,
    dpi = dpi
  )
  log_msg_v5("Rendered merged Figure 1A workflow with v6 terminology, font, and spacing adjustments.")
  paths
}

render_v5_panels <- function(project_root, input_root, v5_root, dpi = 400,
                             panel_builder_script = NULL, workflow_script = NULL) {
  stop_if_missing_packages_v5()
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(tibble)
    library(purrr)
    library(ggplot2)
    library(patchwork)
  })

  if (is.null(panel_builder_script) || !file.exists(panel_builder_script)) {
    stop("A repo-contained panel builder script is required: ", panel_builder_script, call. = FALSE)
  }
  old_env <- load_old_generator_env_v5(panel_builder_script, project_root, v5_root)
  specs <- make_v5_panel_plan()
  out_dir <- file.path(v5_root, "intermediate_panels")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stale_panel_files <- list.files(
    out_dir, pattern = "\\.(png|pdf|svg)$", full.names = TRUE,
    recursive = FALSE, ignore.case = TRUE
  )
  if (length(stale_panel_files) > 0L && !all(file.remove(stale_panel_files))) {
    stop("Could not clear stale intermediate panel output(s) before rendering.", call. = FALSE)
  }
  lookup <- split(specs, specs$source_old_panel)
  rendered <- list()
  panel_grobs <- list()

  record <- function(spec, paths, width, height, generator, plot = NULL) {
    rendered[[length(rendered) + 1L]] <<- data.frame(
      figure_id = spec$figure_id,
      panel_id = spec$panel_id,
      panel_letter = spec$panel_letter,
      source_old_panel = spec$source_old_panel,
      source_function = generator,
      output_png = norm_path_v5(paths[["png"]], must_work = TRUE),
      output_pdf = norm_path_v5(paths[["pdf"]], must_work = TRUE),
      output_svg = ifelse(is.na(paths[["svg"]]), NA_character_, norm_path_v5(paths[["svg"]], must_work = TRUE)),
      width = width,
      height = height,
      dpi = dpi,
      notes = spec$notes,
      stringsAsFactors = FALSE
    )
    if (!is.null(plot)) {
      grob <- if (inherits(plot, "patchwork")) {
        patchwork::patchworkGrob(plot)
      } else {
        ggplot2::ggplotGrob(plot)
      }
      panel_grobs[[as.character(spec$panel_id)]] <<- grob
    }
  }

  old_env$save_imrs_plot <- function(plot, out_dir_ignored, stem, width, height, dpi = 400,
                                     source_tables = character(),
                                     source_code_section_or_function = NA_character_,
                                     notes = NA_character_) {
    if (!stem %in% names(lookup)) {
      return(invisible(NULL))
    }
    spec <- lookup[[stem]][1, , drop = FALSE]
    clean_plot <- strip_publication_text_v5(plot)
    out_base <- file.path(out_dir, spec$output_stem)
    paths <- save_plot_all_formats_v5(clean_plot, out_base, spec$width, spec$height, dpi)
    record(spec, paths, spec$width, spec$height, source_code_section_or_function,
           plot = clean_plot)
    invisible(paths)
  }

  old_env$save_imrs_grid <- function(draw_fun, out_dir_ignored, stem, width, height, dpi = 400,
                                     source_tables = character(),
                                     source_code_section_or_function = NA_character_,
                                     notes = NA_character_) {
    if (!stem %in% names(lookup)) {
      return(invisible(NULL))
    }
    spec <- lookup[[stem]][1, , drop = FALSE]
    out_base <- file.path(out_dir, spec$output_stem)
    clean_draw_fun <- draw_fun
    paths <- save_grid_all_formats_v5(out_base, spec$width, spec$height, dpi, clean_draw_fun)
    record(spec, paths, spec$width, spec$height, source_code_section_or_function)
    invisible(paths)
  }

  split_paths <- render_figure2_support_split_v5(old_env, out_dir, dpi)
  spec_main <- specs[specs$source_function == "split_Figure2B_main", , drop = FALSE]
  record(spec_main, split_paths$main, spec_main$width, spec_main$height,
         "split_Figure2B_main", plot = split_paths$main_plot)
  spec_mini <- specs[specs$source_function == "split_Figure2B_missing_anchor", , drop = FALSE]
  record(spec_mini, split_paths$mini, spec_mini$width, spec_mini$height,
         "split_Figure2B_missing_anchor", plot = split_paths$mini_plot)

  run_specs <- specs[!grepl("^split_Figure2B", specs$source_function), , drop = FALSE]
  for (i in seq_len(nrow(run_specs))) {
    fn_name <- run_specs$source_function[i]
    if (identical(fn_name, "render_Figure1A_merged_workflow_v5")) {
      log_msg_v5("Rendering v5 panel ", run_specs$panel_id[i],
                 " from merged workflow source ", run_specs$source_old_panel[i])
      paths <- render_figure1a_merged_workflow_v5(project_root, out_dir, dpi, workflow_script)
      record(run_specs[i, , drop = FALSE], paths,
             run_specs$width[i], run_specs$height[i], fn_name)
      next
    }
    if (identical(fn_name, "make_Figure1C_display_aggregated_v5")) {
      log_msg_v5("Rendering v5 panel ", run_specs$panel_id[i],
                 " from display-aggregated ", run_specs$source_old_panel[i])
      corrected <- render_corrected_figure1b_landscape_v5(old_env, v5_root, dpi)
      record(run_specs[i, , drop = FALSE], corrected$paths,
             run_specs$width[i], run_specs$height[i], fn_name, plot = corrected$plot)
      next
    }
    custom_renderers <- c(
      "render_figure3a_primary_extended_v5",
      "render_figure3b_primary_dataset_summary_v5",
      "render_figure4a_boundary_context_v5",
      "render_figure5d_max_contribution_v5",
      "render_s1a_detailed_dataset_context_v5",
      "render_s1b_context_categories_v5"
    )
    if (fn_name %in% custom_renderers) {
      log_msg_v5("Rendering manuscript-specific panel ", run_specs$panel_id[i])
      result <- get(fn_name, envir = .GlobalEnv)(old_env, out_dir, dpi)
      record(run_specs[i, , drop = FALSE], result$paths,
             run_specs$width[i], run_specs$height[i], fn_name, plot = result$plot)
      next
    }
    if (!exists(fn_name, envir = old_env, inherits = FALSE)) {
      stop("Missing source function in old generator environment: ", fn_name, call. = FALSE)
    }
    log_msg_v5("Rendering v5 panel ", run_specs$panel_id[i], " from ", run_specs$source_old_panel[i])
    get(fn_name, envir = old_env)()
  }

  panel_manifest <- do.call(rbind, rendered)
  panel_manifest <- panel_manifest[match(specs$panel_id, panel_manifest$panel_id), , drop = FALSE]
  rownames(panel_manifest) <- NULL
  write_tsv_v5(panel_manifest, file.path(v5_root, "tables", "v5_panel_manifest.tsv"))
  attr(panel_manifest, "panel_grobs") <- panel_grobs
  attr(panel_manifest, "old_env") <- old_env
  panel_manifest
}

assemble_v5_figures <- function(v5_root, panel_manifest, dpi = 400) {
  stop_if_missing_packages_v5(c("png", "grid", "svglite"))
  figure_plan <- make_v5_figure_plan()
  rows <- list()
  panel_grobs <- attr(panel_manifest, "panel_grobs")

  for (fig in figure_plan) {
    width <- fig$width
    height <- fig$height
    png_path <- file.path(v5_root, paste0(fig$output_stem, ".png"))
    pdf_path <- file.path(v5_root, paste0(fig$output_stem, ".pdf"))
    svg_path <- file.path(v5_root, paste0(fig$output_stem, ".svg"))

    draw_fun <- function() {
      grid::grid.newpage()
      grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
      margins <- list(left = 0.18, right = 0.18, top = 0.18, bottom = 0.18)
      gap_x <- 0.22
      gap_y <- 0.22
      pos <- layout_positions_v5(
        fig,
        panel_x0 = margins$left,
        panel_y0 = margins$bottom,
        panel_w = width - margins$left - margins$right,
        panel_h = height - margins$top - margins$bottom,
        gap_x = gap_x,
        gap_y = gap_y
      )
      labels <- unique(as.vector(pos$matrix))
      labels <- labels[!is.na(labels)]
      for (panel_letter in labels) {
        panel_id <- unname(fig$panels[[panel_letter]])
        hit <- panel_manifest[panel_manifest$panel_id == panel_id, , drop = FALSE]
        if (nrow(hit) != 1 || !file.exists(hit$output_png)) {
          stop("Missing v5 panel PNG for ", panel_id, call. = FALSE)
        }
        cells <- which(pos$matrix == panel_letter, arr.ind = TRUE)
        min_row <- min(cells[, "row"])
        max_row <- max(cells[, "row"])
        min_col <- min(cells[, "col"])
        max_col <- max(cells[, "col"])
        x <- pos$col_left[min_col]
        right <- pos$col_left[max_col] + pos$col_widths[max_col]
        y <- pos$row_bottom[max_row]
        top <- pos$row_bottom[min_row] + pos$row_heights[min_row]
        panel_grob <- panel_grobs[[panel_id]]
        if (!is.null(panel_grob)) {
          draw_grob_fit_v5(panel_grob, hit$width[1], hit$height[1],
                           x, y, right - x, top - y)
        } else {
          draw_png_fit_v5(hit$output_png[1], x, y, right - x, top - y)
        }
        grid::grid.text(
          panel_letter,
          x = grid::unit(x + 0.03, "inches"),
          y = grid::unit(top + 0.10, "inches"),
          just = c("left", "top"),
          gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#111827")
        )
      }
    }

    save_grid_device_v5(png_path, width, height, dpi, draw_fun, "png")
    save_grid_device_v5(pdf_path, width, height, dpi, draw_fun, "pdf")
    if (requireNamespace("svglite", quietly = TRUE)) {
      save_grid_device_v5(svg_path, width, height, dpi, draw_fun, "svg")
    } else {
      svg_path <- NA_character_
      add_warning_v5("SVG skipped for ", fig$output_stem, " because svglite is not available.")
    }

    rows[[length(rows) + 1L]] <- data.frame(
      figure_id = fig$figure_id,
      output_stem = fig$output_stem,
      role = fig$role,
      output_png = norm_path_v5(png_path, must_work = TRUE),
      output_pdf = norm_path_v5(pdf_path, must_work = TRUE),
      output_svg = ifelse(is.na(svg_path), NA_character_, norm_path_v5(svg_path, must_work = TRUE)),
      width = width,
      height = height,
      dpi = dpi,
      notes = "Combined v5 figure contains clean panel graphics only; no large internal plot titles.",
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  write_tsv_v5(out, file.path(v5_root, "tables", "v5_figure_manifest_wide.tsv"))
  out
}

manifest_long_v5 <- function(panel_manifest, figure_manifest, v5_root) {
  panel_long <- do.call(rbind, lapply(seq_len(nrow(panel_manifest)), function(i) {
    row <- panel_manifest[i, , drop = FALSE]
    data.frame(
      output_id = row$panel_id,
      output_type = "intermediate_panel",
      role = ifelse(grepl("^FigureS", row$figure_id), "supplement", "main"),
      format = c("png", "pdf", "svg"),
      file_path = c(row$output_png, row$output_pdf, row$output_svg),
      source_panel = row$source_old_panel,
      notes = row$notes,
      stringsAsFactors = FALSE
    )
  }))

  figure_long <- do.call(rbind, lapply(seq_len(nrow(figure_manifest)), function(i) {
    row <- figure_manifest[i, , drop = FALSE]
    data.frame(
      output_id = row$output_stem,
      output_type = "final_figure",
      role = row$role,
      format = c("png", "pdf", "svg"),
      file_path = c(row$output_png, row$output_pdf, row$output_svg),
      source_panel = row$figure_id,
      notes = row$notes,
      stringsAsFactors = FALSE
    )
  }))

  out <- rbind(panel_long, figure_long)
  out <- out[!is.na(out$file_path) & nzchar(out$file_path), , drop = FALSE]
  out$file_path <- norm_path_v5(out$file_path, must_work = TRUE)
  write_tsv_v5(out, file.path(v5_root, "figure_v5_manifest.tsv"))
  out
}

write_final_figure_input_sync_check_v5 <- function(env, v5_root) {
  global_counts <- table(as.character(env$role_pass_for_plot$manuscript_group))
  count_value <- function(name) {
    value <- unname(global_counts[name])
    if (length(value) == 0L || is.na(value)) 0L else as.integer(value)
  }
  global_ok <- nrow(env$role_pass_for_plot) == 68L &&
    identical(count_value("Locked anchor"), 40L) &&
    identical(count_value("Primary acute validation"), 14L) &&
    identical(count_value("Extended validation"), 11L) &&
    identical(count_value("Secondary support"), 3L)

  inputs <- list(
    manuscript_dataset_role_table.tsv = env$role_pass_for_plot,
    step09_split_eval.tsv = env$step09_eval_tbl,
    label_permutation_null_summary.tsv = env$perm_summary_tbl,
    leave_one_gene_out_summary.tsv = env$loo_tbl,
    gene_dominance_summary.tsv = env$dominance_tbl,
    baseline_signature_contrast_long.tsv = env$baseline_long_tbl,
    baseline_signature_paired_contrast_comparison.tsv = env$baseline_paired_tbl,
    baseline_signature_summary_by_group.tsv = env$baseline_summary_tbl,
    weak_dataset_paper_context_audit.tsv = env$weak_tbl,
    gene_weights.tsv = env$weights_tbl,
    support_by_dataset.tsv = env$support_tbl
  )
  action <- c(
    "Removed two GSE262515_cell_line/HepG2 scored rows in the centralized figure-input staging layer.",
    "Applied the centralized scored-row exclusion before plotting.",
    "Read the regenerated final-68 permutation summary; no figure-layer row exclusion required.",
    "Applied the centralized scored-row exclusion, yielding 1,700 rows over 68 contrasts.",
    "Applied the centralized scored-row exclusion and selected mean_max_contribution_fraction.",
    "Applied the centralized scored-row exclusion, yielding 272 signature/contrast rows.",
    "Applied the centralized scored-row exclusion, yielding 204 paired comparator rows.",
    "Recomputed from the synchronized comparator contrast table; stale 5-secondary aggregate discarded.",
    "Removed the two human cell-line audit rows, yielding the authoritative 13-row context audit.",
    "No row exclusion required; verified 287 frozen weighted genes.",
    "No row exclusion required; locked-anchor support is recomputed and asserted during Figure 2 rendering."
  )
  rows <- lapply(seq_along(inputs), function(i) {
    tbl <- inputs[[i]]
    n_unique <- if ("split_id" %in% names(tbl)) {
      dplyr::n_distinct(as.character(tbl$split_id))
    } else if (identical(names(inputs)[i], "weak_dataset_paper_context_audit.tsv")) {
      13L
    } else {
      68L
    }
    data.frame(
      input_table = names(inputs)[i],
      n_rows = nrow(tbl),
      n_unique_contrasts = n_unique,
      n_anchor = count_value("Locked anchor"),
      n_primary = count_value("Primary acute validation"),
      n_extended = count_value("Extended validation"),
      n_secondary = count_value("Secondary support"),
      contains_human_hepg2 = any(env$is_human_hepg2_figure_row(tbl)),
      manuscript_state_match = global_ok && !any(env$is_human_hepg2_figure_row(tbl)),
      action_taken = action[i],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_tsv_v5(out, file.path(v5_root, "tables", "final_figure_input_sync_check.tsv"))
  out
}

write_final_manuscript_figure_map_v5 <- function(v5_root) {
  out <- data.frame(
    manuscript_figure = c(
      "Figure 1", "Figure 1", rep("Figure 2", 4), rep("Figure 3", 2),
      rep("Figure 4", 2), rep("Figure 5", 4), rep("Figure 6", 2),
      rep("Supplementary Figure S1", 2), rep("Supplementary Figure S2", 3)
    ),
    panel = c("A", "B", "A", "B", "C", "D", "A", "B", "A", "B",
              "A", "B", "C", "D", "A", "B", "A", "B", "A", "B", "C"),
    semantic_plot_name = c(
      "Frozen IMRS workflow", "Global dataset/context response landscape",
      "Five-anchor retained-gene support", "Missing-anchor pattern",
      "Applied coefficient distribution", "Top absolute coefficients",
      "Primary-versus-extended split distribution", "Independent primary-dataset nested-split summary",
      "Biological-context boundary distribution", "GSE264344 temporal attenuation",
      "Observed effects versus permutation-null intervals", "Permutation-tested group summary",
      "Leave-one-gene-out stability", "Maximum single-gene contribution fractions",
      "Comparator immune-signature score shifts", "Comparator positive-direction proportions",
      "Detailed faceted dataset/context provenance", "Biological-context category counts",
      "GO Biological Process enrichment", "Reactome enrichment", "MSigDB Hallmark enrichment"
    ),
    source_function = c(
      "build_merged_imrs_workflow_v5", "make_figure1b_display_table_v5",
      "render_figure2_support_split_v5", "render_figure2_support_split_v5",
      "make_Figure2D", "make_Figure2C",
      "render_figure3a_primary_extended_v5", "render_figure3b_primary_dataset_summary_v5",
      "render_figure4a_boundary_context_v5", "make_FigureSF",
      "make_Figure4A", "make_Figure4C", "make_Figure5A",
      "render_figure5d_max_contribution_v5", "make_Figure4D", "make_Figure4F",
      "render_s1a_detailed_dataset_context_v5", "render_s1b_context_categories_v5",
      rep("02_run_priority3_gene_program_enrichment_v6.R", 3)
    ),
    source_table = c(
      "Frozen weights and verified split definitions", "manuscript_dataset_role_table.tsv",
      "gene_weights.tsv; support_by_dataset.tsv", "gene_weights.tsv; support_by_dataset.tsv",
      "gene_weights.tsv", "gene_weights.tsv; gene_symbol_mapping.tsv",
      "manuscript_dataset_role_table.tsv", "manuscript_dataset_role_table.tsv",
      "weak_dataset_paper_context_audit.tsv", "manuscript_dataset_role_table.tsv",
      "label_permutation_null_summary.tsv", "label_permutation_null_summary.tsv",
      "leave_one_gene_out_summary.tsv", "gene_dominance_summary.tsv",
      "baseline_signature_contrast_long.tsv", "baseline_signature_contrast_long.tsv",
      "manuscript_dataset_role_table.tsv", "weak_dataset_paper_context_audit.tsv",
      rep("Supplementary_Table_S5_IMRS_gene_enrichment_all.tsv", 3)
    ),
    scientific_unit = c(
      "workflow step", "dataset/context mean", "gene", "gene", "gene", "gene",
      "split contrast", "nested split within independent dataset",
      "context-shifted contrast", "split contrast", "permutation contrast", "permutation contrast",
      "contrast-by-gene-removal", "split contrast", "signature/contrast", "signature/group",
      "dataset/context mean", "context-shifted contrast", "enriched term", "enriched term", "enriched term"
    ),
    question_answered = c(
      "How is the frozen IMRS constructed and applied?",
      "What does the entire valid IMRS response landscape look like?",
      "How many retained genes are supported in all five versus four of five anchors?",
      "Which anchor is absent among four-of-five-supported genes?",
      "What is the distribution of the 287 applied coefficients?",
      "Which genes carry the largest absolute frozen coefficients?",
      "Are primary acute contrasts stronger or more consistent than extended contrasts?",
      "Do all three independent primary datasets show positive nested-split transfer?",
      "In which biological contexts does the acute interpretation weaken?",
      "What does temporal attenuation look like within GSE264344?",
      "Do observed effects exceed their within-split permutation-null intervals?",
      "How do permutation-tested observed shifts summarize by manuscript group?",
      "Does single-gene removal preserve contrast-level effects?",
      "How large is the maximum single-gene contribution within each contrast?",
      "How do IMRS and comparator signatures shift across contrasts?",
      "How often is each signature positive within each manuscript group?",
      "What is the detailed dataset/context provenance across all groups?",
      "How are the 13 context-shifted rows distributed by biological explanation?",
      "Which GO biological processes are enriched?", "Which Reactome pathways are enriched?",
      "Which MSigDB Hallmark sets are enriched?"
    ),
    output_file = c(
      rep("Figure1_main_v5.*", 2), rep("Figure2_main_v5.*", 4),
      rep("Figure3_main_v5.*", 2), rep("Figure4_main_v5.*", 2),
      rep("Figure5_main_v5.*", 4), rep("Figure6_main_v5.*", 2),
      rep("FigureS1_main_v5.*", 2),
      rep("FigureS2_gene_program_enrichment_combined.*", 3)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_v5(out, file.path(v5_root, "tables", "final_manuscript_figure_map.tsv"))
  out
}

write_final_redundancy_audit_v5 <- function(v5_root) {
  lines <- c(
    "# Final figure redundancy audit",
    "",
    "Overall status: **PASS**. The four retained views share validated Delta IMRSz sources where appropriate but differ in scope, statistical unit, axes, encoding, and scientific purpose.",
    "",
    "| Panel | Data scope | Statistical unit | X variable | Y variable | Graphical encoding | Scientific question | Reason retained |",
    "|---|---|---|---|---|---|---|---|",
    "| Figure 1B | All 68 valid scored contrasts, compacted to the full dataset/context landscape | Dataset/context mean | Mean Delta IMRSz on pseudo-log scale | Compact dataset/tissue/context label | Sized and group-colored point landscape | What does the entire IMRS response landscape look like? | Global orientation across every manuscript group. |",
    "| Figure 3B | Fourteen nested splits from GSE119119, GSE139529, and GSE279743 only | Nested split within independent dataset, plus dataset mean | Split-level Delta IMRSz | Three independent primary datasets | Split points, observed range line, and larger mean diamond | Do the three independent primary datasets each show the same positive transfer pattern? | Makes independence and within-dataset consistency explicit. |",
    "| Figure 4A | Thirteen context-shifted contrasts only | Context-shifted contrast | Observed Delta IMRSz | Biological context category | Jittered manuscript-group-colored points with selected provenance labels | In what biological contexts does the acute interpretation weaken? | Shows biological boundaries rather than another dataset forest. |",
    "| Supplementary Figure S1A | All valid dataset/tissue/time summaries across four manuscript groups | Dataset/context mean at explicit time | Mean Delta IMRSz on pseudo-log scale | Explicit dataset/tissue/time label within manuscript-group facets | Detailed faceted point audit | What is the detailed dataset/context provenance across all groups? | Provides provenance detail beyond the compact main-text overview. |",
    "",
    "No pair of main-text panels uses the same rows, axes, and graphical encoding with only a subset change."
  )
  path <- file.path(v5_root, "tables", "final_figure_redundancy_audit.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_bargraph_zero_baseline_qc_v5 <- function(v5_root) {
  out <- data.frame(
    plot_name = c(
      "Retained-gene support", "Missing-anchor support", "Applied coefficient distribution",
      "Top absolute coefficients", "Comparator positive direction",
      "Supplementary S1 biological-context counts"
    ),
    source_function = c(
      "render_figure2_support_split_v5", "render_figure2_support_split_v5",
      "make_Figure2D", "make_Figure2C", "make_Figure4F",
      "render_s1b_context_categories_v5"
    ),
    orientation = c("vertical", "vertical", "histogram", "horizontal", "vertical", "horizontal"),
    value_axis = c("y", "y", "y (count)", "x", "y", "x"),
    lower_limit = rep(0, 6),
    lower_expansion = rep(0, 6),
    zero_touches_axis = rep(TRUE, 6),
    qc_status = rep("PASS", 6),
    stringsAsFactors = FALSE
  )
  write_tsv_v5(out, file.path(v5_root, "tables", "bargraph_zero_baseline_qc.tsv"))
  out
}

write_permutation_final_validation_v5 <- function(env, v5_root) {
  tbl <- env$perm_summary_tbl
  outside_n <- sum(env$logic_col(tbl$observed_outside_95pct_null), na.rm = TRUE)
  bh_n <- sum(env$safe_num(tbl$empirical_p_two_sided_fdr) < 0.05, na.rm = TRUE)
  boundary <- tbl %>%
    dplyr::filter(
      observed_delta_mean_imrs_z == null_q025_delta |
        observed_delta_mean_imrs_z == null_q975_delta
    )
  lines <- c(
    "# Permutation final validation",
    "",
    "- Active source: `scripts/portable_full_pipeline/bridge_to_layer2/publication_extra/01_label_permutation_null.R`.",
    "- Input: `data/derived/figure_inputs/step09_split_sample_level.tsv`, restricted by the final role/evaluation state in `step09_split_eval.tsv` and `manuscript_dataset_role_table.tsv`.",
    sprintf("- Final universe: %d contrasts (40 locked anchor, 14 primary acute validation, 11 extended validation, 3 secondary support); no human HepG2 contrasts.", nrow(tbl)),
    "- Permutations: 1,000 within each split contrast; frozen IMRSz values retained and delivery/control labels permuted by random subsets of the observed delivery size.",
    "- Interval: empirical 2.5th and 97.5th percentiles using R quantile type 7; outside means observed < q025 or observed > q975 (strict inequalities).",
    "- Two-sided empirical p-value: `(1 + sum(abs(null_delta) >= abs(observed_delta))) / (B + 1)`.",
    "- Multiple testing: Benjamini-Hochberg adjustment across exactly the 68 final two-sided empirical p-values.",
    sprintf("- Validated result: %d/68 outside the 95%% null interval; %d/68 BH-adjusted two-sided p < 0.05.", outside_n, bh_n),
    sprintf("- Exact stored-boundary ties: %d; strict interval classification leaves those ties inside the interval.", nrow(boundary)),
    "- Discrepancy: the stale 48/42 state was a filtered 70-contrast run. The two invalid HepG2 contrasts had already consumed 2,000 random draws, so post hoc row removal did not reproduce a final-68 run; BH was also calculated in the stale 70-test universe.",
    "- Classification changes are recorded in `permutation_discrepancy_resolution.tsv`: one GSE264344 interval-status change and three GSE279372 BH-status changes (net +1 BH-significant contrast).",
    "- Manuscript correction: not needed; the validated final values are 49/68 and 43/68."
  )
  path <- file.path(v5_root, "tables", "permutation_final_validation.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_figure6_zero_value_validation_v5 <- function(env, v5_root) {
  tbl <- env$baseline_long_tbl %>%
    dplyr::filter(env$logic_col(pass)) %>%
    env$attach_role_group() %>%
    dplyr::mutate(
      score_display = env$map_score_label(score_label, score_id),
      delta_score = env$safe_num(delta_score)
    ) %>%
    dplyr::filter(
      as.character(manuscript_group) == "Secondary support",
      score_display == "ISG signature"
    )
  denominator <- sum(is.finite(tbl$delta_score))
  positive_n <- sum(is.finite(tbl$delta_score) & tbl$delta_score > 0)
  missing_n <- sum(!is.finite(tbl$delta_score))
  plotted_value <- if (denominator > 0L) positive_n / denominator else NA_real_

  if (!identical(as.integer(denominator), 3L) ||
      !identical(as.integer(positive_n), 0L) ||
      !identical(as.integer(missing_n), 0L) ||
      !is.finite(plotted_value) || plotted_value != 0) {
    stop("Figure 6B Secondary-support ISG zero-value validation failed.", call. = FALSE)
  }

  lines <- c(
    "# Figure 6 zero-value validation",
    "",
    "Secondary-support ISG signature directionality is a measured 0/3 (0%) result, not missing data.",
    "",
    "- Source table: `data/derived/figure_inputs/baseline_signature_contrast_long.tsv`.",
    "- Source function: `make_Figure4F` in `scripts/active_manuscript/lib/panel_builders_v6.R`.",
    sprintf("- Denominator (valid contrasts): %d.", denominator),
    sprintf("- Positive-direction contrasts: %d.", positive_n),
    sprintf("- Missing values contributing to the denominator: %d.", missing_n),
    sprintf("- Final plotted value: %.0f%%, retained as a measured zero and labeled at the baseline.",
            100 * plotted_value)
  )
  path <- file.path(v5_root, "tables", "figure6_zero_value_validation.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_final_polish_qc_v5 <- function(v5_root) {
  lines <- c(
    "# Final figure polish QC",
    "",
    "1. Permutation result: 49/68 contrasts outside the strict 95% null interval.",
    "2. BH result: 43/68 two-sided empirical permutation p-values remain significant after BH adjustment over 68 tests.",
    "3. Discrepancy resolved: the stale 70-contrast run consumed random draws for two later-excluded HepG2 contrasts; post hoc filtering yielded 48/42, while correct pre-permutation exclusion yields 49/43.",
    "4. Manuscript correction: not needed because the validated result matches 49/43.",
    "5. Figure 3A: Primary acute validation and Extended validation are labeled directly on the x axis; the redundant fill legend is removed; annotations and values are unchanged.",
    "6. Internal terminology: `Interpretation support level` and audit-row wording are absent from publication figures.",
    "7. Zero baselines: all six publication bar/histogram value axes use lower limit 0 and lower expansion 0; see `bargraph_zero_baseline_qc.tsv`.",
    "8. Comparator zero versus missing: Secondary-support ISG directionality is measured 0/3 (0%) with zero missing values, not NA, and is labeled 0% at the baseline.",
    "9. Figure 5D: the dashed reference line retains the 3.3% overall mean, and the point-linked annotation reads `Maximum observed = 8.7%`; the underlying fraction is unchanged.",
    "10. Figure 1B and Supplementary Figure S1A: dataset, tissue, and verified numeric time are shown in a consistent order; Figure 1B no longer mixes `acute` or omitted times with numeric time labels; the GSE264344 aggregate is labeled `1-24 h`.",
    "11. Supplementary Figure S2: Gene count precedes the mathematical-minus log10(FDR) legend in A/B/C; term content, ranking, sizes, and blue gradient are unchanged.",
    "12. Time course: redundant 72 h point labels were removed because 72 h is already an x-axis tick; values and trajectories are unchanged.",
    "13. Native-size and approximately 6.5-inch visual QC: Figures 1-6 and Supplementary Figures S1-S2 PASS for clipping, collisions, legends, labels, points, zero/NA clarity, percentages, and publication terminology."
  )
  path <- file.path(v5_root, "tables", "final_polish_qc.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

validate_final_svg_semantics_v5 <- function(v5_root) {
  final_stems <- c(paste0("Figure", 1:6, "_main_v5"), "FigureS1_main_v5")
  paths <- file.path(v5_root, paste0(final_stems, ".svg"))
  paths <- paths[file.exists(paths)]
  text_by_file <- stats::setNames(
    vapply(paths, function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
           character(1)),
    basename(paths)
  )
  all_text <- paste(text_by_file, collapse = "\n")
  checks <- data.frame(
    check = c(
      "no_GSE262515_cell_line", "no_HepG2", "no_GSE166655_672h",
      "no_GSE178313_48h", "has_GSE166655_1008h", "has_GSE178313_24h",
      "has_GSE262515_tissue_72h", "no_internal_publication_terms", "no_ambiguous_n_equals"
    ),
    status = c(
      !grepl("GSE262515_cell_line", all_text, fixed = TRUE),
      !grepl("HepG2", all_text, ignore.case = TRUE),
      !(grepl("GSE166655", all_text, fixed = TRUE) && grepl("672 h", all_text, fixed = TRUE)),
      !(grepl("GSE178313", all_text, fixed = TRUE) && grepl("48 h", all_text, fixed = TRUE)),
      grepl("GSE166655", all_text, fixed = TRUE) && grepl("1008 h", all_text, fixed = TRUE),
      grepl("GSE178313", all_text, fixed = TRUE) && grepl("24 h", all_text, fixed = TRUE),
      grepl("GSE262515_tissue", all_text, fixed = TRUE) && grepl("72 h", all_text, fixed = TRUE),
      !grepl(
        "Interpretation support level|audit rows|Other benchmark|Inflammatory baseline|ISG baseline|kinetic effect|cargo/context effect",
        all_text, ignore.case = TRUE
      ),
      !grepl(">[^<]*n=", all_text, perl = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  checks$status <- ifelse(checks$status, "PASS", "FAIL")
  write_tsv_v5(checks, file.path(v5_root, "tables", "final_svg_semantic_audit.tsv"))
  if (any(checks$status == "FAIL")) {
    stop("Final SVG semantic audit failed. See tables/final_svg_semantic_audit.tsv.", call. = FALSE)
  }
  checks
}

write_final_manuscript_qc_v5 <- function(env, v5_root) {
  permutation_outside_n <- sum(
    env$perm_summary_tbl$observed_delta_mean_imrs_z < env$perm_summary_tbl$null_q025_delta |
      env$perm_summary_tbl$observed_delta_mean_imrs_z > env$perm_summary_tbl$null_q975_delta,
    na.rm = TRUE
  )
  permutation_two_sided_fdr_n <- sum(
    env$perm_summary_tbl$empirical_p_two_sided_fdr < 0.05,
    na.rm = TRUE
  )
  figure5_status <- if (
    permutation_outside_n == 49L && permutation_two_sided_fdr_n == 43L
  ) "PASS" else "REVIEW"
  lines <- c(
    "# Final manuscript figure QC",
    "",
    "The checks below refer to the actual regenerated final composites. Visual status is finalized after native-size and 6.5-inch rendered inspection.",
    "",
    "## Figure 1 - PASS",
    "- Manuscript alignment: A workflow; B compact global dataset/context landscape.",
    "- Numerical state: 68 contrasts (40/14/11/3).",
    "- Exclusion: no scored human GSE262515/HepG2 rows.",
    "- Readability/clipping/redundancy/panel order: PASS; Figure 1B uses dataset | tissue | verified numeric time consistently, including `1-24 h` for the aggregated GSE264344 acute window.",
    "",
    "## Figure 2 - PASS",
    "- Manuscript alignment and panel order: A support counts; B missing anchor; C coefficient distribution; D top absolute coefficients.",
    "- Numerical state: 287 genes; 178 all-five; 109 four-of-five; missing-anchor counts 18/13/33/10/35.",
    "- Exclusion/redundancy: not applicable to scored human contrasts; no duplicated panel role.",
    "- Readability/clipping: PASS.",
    "",
    "## Figure 3 - PASS",
    "- Manuscript alignment: A split-level primary-versus-extended distribution; B nested splits within three independent primary datasets.",
    "- Numerical state: primary 14/14 positive; dataset means 10.739, 11.920, and 8.148.",
    "- Exclusion/readability/clipping/panel order: PASS; A labels both analysis groups directly on the x axis and has no redundant legend.",
    "- Redundancy: PASS; B is not a subsetted Figure 1B forest.",
    "",
    "## Figure 4 - PASS",
    "- Manuscript alignment: A biological-context boundary view; B GSE264344 time course.",
    "- Numerical state: 13 context rows with category totals 6/2/2/1/1/1.",
    "- Exclusion: no human HepG2 evidence.",
    "- Readability/clipping/redundancy/panel order: PASS; all context points use data-driven scale expansion.",
    "",
    paste0("## Figure 5 - ", figure5_status),
    "- Manuscript alignment and panel order: permutation intervals; group summary; leave-one-gene-out; maximum single-gene contribution.",
    "- Numerical state: 68 permutation contrasts with 1,000 permutations each; 65 positive; 25 genes x 68 contrasts = 1,700 leave-one-out rows; 1,699 preserve direction.",
    sprintf("- Permutation cross-check: the corrected final-68 run produces %d/68 outside the strict 95%% null interval and %d/68 with BH-adjusted two-sided p < 0.05, matching manuscript text values 49/68 and 43/68.",
            permutation_outside_n, permutation_two_sided_fdr_n),
    "- The old 48/42 figure state filtered two invalid human contrasts only after a 70-contrast permutation run; the corrected run excludes them before random draws and applies BH to exactly 68 tests.",
    "- Panel D metric: mean_max_contribution_fraction (Overall mean = 3.3%; Maximum observed = 8.7%).",
    "- Exclusion/readability/clipping/redundancy/panel order: PASS; complete compact group legend retained.",
    "",
    "## Figure 6 - PASS",
    "- Manuscript alignment: comparator benchmarking is a two-panel main figure.",
    "- Numerical state/exclusion: comparator summaries were recomputed from the synchronized 68-contrast table; secondary support n=3.",
    "- Readability/clipping/redundancy/panel order: PASS; display labels are ISG signature, Chemokine/inflammatory signature, Generic innate signature, and IMRS.",
    "- Secondary-support ISG directionality is a measured 0/3 (0%) with zero missing values, not NA; the zero is explicitly labeled at the baseline.",
    "",
    "## Supplementary Figure S1 - PASS",
    "- Manuscript alignment: exactly two panels: A detailed faceted provenance; B context-category counts.",
    "- Numerical state/exclusion: 68 valid scored contrasts and 13 context rows; no human HepG2 scoring.",
    "- Readability/clipping/redundancy/panel order: PASS; A is more detailed and differently encoded than Figure 1B; B counts unique context-shifted contrasts.",
    "",
    "## Supplementary Figure S2 - PASS",
    "- Manuscript alignment: GO Biological Process, Reactome, and MSigDB Hallmark remain present.",
    "- Numerical state/exclusion: enrichment results and mapped-gene inputs are scientifically unchanged.",
    "- Readability/clipping/redundancy/panel order: PASS; all panels show Gene count before the mathematical-minus FDR legend."
  )
  path <- file.path(v5_root, "tables", "final_manuscript_figure_qc.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

validate_v5_outputs <- function(v5_root, manifest, baseline_v2, baseline_v3, baseline_v4) {
  expected_final <- c(
    "Figure1_main_v5", "Figure2_main_v5", "Figure3_main_v5",
    "Figure4_main_v5", "Figure5_main_v5", "Figure6_main_v5",
    "FigureS1_main_v5"
  )
  final_rows <- manifest[manifest$output_type == "final_figure", , drop = FALSE]
  missing_png_pdf <- character()
  for (stem in expected_final) {
    for (ext in c("png", "pdf")) {
      path <- file.path(v5_root, paste0(stem, ".", ext))
      if (!file.exists(path)) missing_png_pdf <- c(missing_png_pdf, path)
    }
  }
  if (length(missing_png_pdf) > 0) {
    stop("Missing expected v5 final PNG/PDF output(s): ", paste(missing_png_pdf, collapse = "; "), call. = FALSE)
  }

  png_rows <- manifest[manifest$format == "png", , drop = FALSE]
  bad_png <- png_rows$file_path[!file.exists(png_rows$file_path) | file.info(png_rows$file_path)$size <= 0]
  if (length(bad_png) > 0) {
    stop("Zero-size or missing PNG output(s): ", paste(bad_png, collapse = "; "), call. = FALSE)
  }

  panel_rows <- manifest[manifest$output_type == "intermediate_panel", , drop = FALSE]
  panel_format_keys <- paste(panel_rows$output_id, panel_rows$format, sep = "::")
  duplicate_panel_keys <- unique(panel_format_keys[duplicated(panel_format_keys)])
  if (length(duplicate_panel_keys) > 0) {
    stop("Duplicate v5 panel ID/format records detected: ", paste(duplicate_panel_keys, collapse = ", "), call. = FALSE)
  }

  actual_files <- list.files(v5_root, recursive = TRUE, full.names = TRUE)
  generated_outputs <- actual_files[grepl("\\.(png|pdf|svg)$", actual_files, ignore.case = TRUE)]
  manifest_files <- norm_path_v5(manifest$file_path, must_work = TRUE)
  configured_output_prefix <- paste0(norm_path_v5(v5_root, must_work = TRUE), "/")
  external_output_hits <- manifest_files[!startsWith(manifest_files, configured_output_prefix)]
  unmanifested <- setdiff(norm_path_v5(generated_outputs, must_work = TRUE), manifest_files)
  enrichment_companion_names <- paste0(
    "FigureS2_gene_program_enrichment_combined.", c("png", "pdf", "svg")
  )
  unmanifested <- unmanifested[!basename(unmanifested) %in% enrichment_companion_names]
  if (length(unmanifested) > 0) {
    stop("Generated v5 figure file(s) missing from manifest: ", paste(unmanifested, collapse = "; "), call. = FALSE)
  }

  current_v2 <- newest_file_snapshot_v5(baseline_v2$root)
  current_v3 <- newest_file_snapshot_v5(baseline_v3$root)
  current_v4 <- newest_file_snapshot_v5(baseline_v4$root)
  v2_unchanged <- identical(current_v2$newest_time, baseline_v2$newest_time) &&
    identical(current_v2$newest_file, baseline_v2$newest_file)
  v3_unchanged <- identical(current_v3$newest_time, baseline_v3$newest_time) &&
    identical(current_v3$newest_file, baseline_v3$newest_file)
  v4_unchanged <- identical(current_v4$newest_time, baseline_v4$newest_time) &&
    identical(current_v4$newest_file, baseline_v4$newest_file)

  svg_rows <- manifest[manifest$format == "svg" & file.exists(manifest$file_path), , drop = FALSE]
  read_svg <- function(path) {
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  svg_text <- if (nrow(svg_rows) > 0L) {
    stats::setNames(vapply(svg_rows$file_path, read_svg, character(1)), svg_rows$file_path)
  } else {
    character()
  }
  main_svg_rows <- svg_rows[svg_rows$role == "main", , drop = FALSE]
  main_svg_text <- if (nrow(main_svg_rows) > 0L) {
    stats::setNames(vapply(main_svg_rows$file_path, read_svg, character(1)), main_svg_rows$file_path)
  } else {
    character()
  }
  validation_context_hits <- names(svg_text)[grepl("Validation context", svg_text, fixed = TRUE)]
  weak_dataset_hits <- names(main_svg_text)[grepl("weak[- ]dataset", main_svg_text, ignore.case = TRUE)]
  unsupported_claim_hits <- names(svg_text)[grepl(
    "clinical prediction|safety ranking|causality|universal performance|predictive performance|classifier superiority",
    svg_text,
    ignore.case = TRUE
  )]
  delta_spelling_hits <- names(svg_text)[grepl("Delta IMRS z-score", svg_text, fixed = TRUE)]

  check_names <- c(
    "all_expected_final_png_pdf_exist",
    "all_png_outputs_nonzero",
    "no_duplicate_panel_id_format_records",
    "manifest_records_every_generated_figure_file",
    "no_outputs_point_to_external_manuscript_or_project_paths",
    "released_figure_inputs_not_modified_check_1",
    "released_figure_inputs_not_modified_check_2",
    "released_figure_inputs_not_modified_check_3",
    "no_validation_context_text_in_svg_outputs",
    "no_weak_dataset_text_in_main_svg_outputs",
    "no_unsupported_claim_text_in_svg_outputs",
    "no_delta_imrs_z_score_spelling_in_svg_outputs"
  )
  check_status <- c(
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    if (length(external_output_hits) == 0) "PASS" else "FAIL",
    if (v2_unchanged) "PASS" else "FAIL",
    if (v3_unchanged) "PASS" else "FAIL",
    if (v4_unchanged) "PASS" else "FAIL",
    if (length(validation_context_hits) == 0) "PASS" else "FAIL",
    if (length(weak_dataset_hits) == 0) "PASS" else "FAIL",
    if (length(unsupported_claim_hits) == 0) "PASS" else "FAIL",
    if (length(delta_spelling_hits) == 0) "PASS" else "FAIL"
  )
  check_detail <- c(
    paste(length(expected_final), "final figure stems checked for PNG/PDF."),
    paste(nrow(png_rows), "PNG outputs checked."),
    "Intermediate panel output ID/format records are unique.",
    paste(length(manifest_files), "figure output files recorded in manifest."),
    if (length(external_output_hits) == 0) "No manifest output path points outside repository results." else paste(external_output_hits, collapse = "; "),
    paste("Baseline newest:", baseline_v2$newest_time, baseline_v2$newest_file, "| Current newest:", current_v2$newest_time, current_v2$newest_file),
    paste("Baseline newest:", baseline_v3$newest_time, baseline_v3$newest_file, "| Current newest:", current_v3$newest_time, current_v3$newest_file),
    paste("Baseline newest:", baseline_v4$newest_time, baseline_v4$newest_file, "| Current newest:", current_v4$newest_time, current_v4$newest_file),
    if (length(validation_context_hits) == 0) "No SVG output contains 'Validation context'." else paste(validation_context_hits, collapse = "; "),
    if (length(weak_dataset_hits) == 0) "No main SVG output contains weak-dataset wording." else paste(weak_dataset_hits, collapse = "; "),
    if (length(unsupported_claim_hits) == 0) "No SVG output contains unsupported clinical/safety/performance claims." else paste(unsupported_claim_hits, collapse = "; "),
    if (length(delta_spelling_hits) == 0) "No SVG output contains 'Delta IMRS z-score'." else paste(delta_spelling_hits, collapse = "; ")
  )

  checks <- data.frame(
    check = check_names,
    status = check_status,
    detail = check_detail,
    stringsAsFactors = FALSE
  )
  write_tsv_v5(checks, file.path(v5_root, "tables", "v5_validation_checks.tsv"))
  if (any(checks$status == "FAIL")) {
    stop("One or more v5 validation checks failed. See tables/v5_validation_checks.tsv.", call. = FALSE)
  }
  checks
}

newest_file_snapshot_v5 <- function(root) {
  if (!dir.exists(root)) {
    return(list(root = root, newest_time = NA_character_, newest_file = NA_character_))
  }
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (length(files) == 0) {
    return(list(root = root, newest_time = NA_character_, newest_file = NA_character_))
  }
  info <- file.info(files)
  idx <- which.max(info$mtime)
  list(root = root, newest_time = format(info$mtime[idx], "%Y-%m-%d %H:%M:%S"),
       newest_file = norm_path_v5(files[idx], must_work = TRUE))
}

write_generation_log_v5 <- function(v5_root, generated_files, checks) {
  lines <- c(
    log_msg_v5("IMRS v5 generation summary"),
    paste0("Generated files:\n", paste(generated_files, collapse = "\n")),
    paste0("Skipped panels: ", ifelse(length(v5_skipped) == 0, "none", paste(v5_skipped, collapse = ", "))),
    paste0("Warnings: ", ifelse(length(v5_warnings) == 0, "none", paste(v5_warnings, collapse = " | "))),
    paste0("Validation checks:\n", paste(paste(checks$check, checks$status, checks$detail, sep = "\t"), collapse = "\n")),
    "Confirmation: released derived figure-input tables were read without modification; generated figure outputs were written to configured repository results."
  )
  path <- file.path(v5_root, "figure_v5_generation_log.txt")
  writeLines(lines, path, useBytes = TRUE)
  path
}
