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
  "Comparator signature" = "Immune-response signature",
  "Benchmark score" = "Immune-response signature",
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
        nrow = 1,
        order = 1,
        override.aes = list(color = "black")
      ),
      color = guide_legend(
        title = "Manuscript analysis group",
        nrow = 2,
        byrow = TRUE,
        order = 2,
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
    theme_imrs_publication(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9.8),
      axis.text.x = element_text(size = 9.2),
      axis.title = element_text(size = 10.2),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.margin = margin(t = 4, r = 12, b = 4, l = 12),
      plot.margin = margin(10, 28, 16, 12)
    )
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1C_simplified_by_dataset", 8.2, 5.5, dpi = 400,
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
        nrow = 1,
        order = 1,
        override.aes = list(color = "black")
      ),
      color = guide_legend(
        title = "Manuscript analysis group",
        nrow = 2,
        byrow = TRUE,
        order = 2,
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
    theme_imrs_publication(base_size = 9.6) +
    theme(
      axis.text.y = element_text(size = 9.6),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8.8),
      legend.title = element_text(size = 9.3),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.margin = margin(t = 4, r = 12, b = 4, l = 12),
      plot.margin = margin(10, 28, 16, 12)
    )
  save_imrs_plot(p, folder_path("FigureS1_weak_late_context_summary"),
                 "FigureS1D_weak_zoom_forest", 8.8, 4.9, dpi = 400,
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

make_v5_panel_plan <- function() {
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
      7.2, 5.8, 3.6, 7.0,
      7.2, 9.5, 10.0,
      9.8, 7.5, 9.0,
      8.4, 8.2, 8.8, 8.6,
      9.0, 8.8
    ),
    height = c(
      6.8, 6.8,
      5.2, 5.2, 5.2, 4.6,
      4.8, 5.5, 6.8,
      4.9, 4.6, 5.5,
      5.1, 5.2, 5.2, 5.1,
      6.4, 5.2
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

make_v5_figure_plan <- function() {
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
      col_widths = c(7.2, 5.8, 3.6),
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

  main_y_max <- max(plot_tbl$n_retained_genes, na.rm = TRUE)
  main_plot <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = support_category, y = n_retained_genes)) +
    ggplot2::geom_col(fill = "#4E79A7", color = "#374151", linewidth = 0.45, width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = n_retained_genes), vjust = -0.35, size = 4.2, color = "#111827") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.16)), breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::coord_cartesian(ylim = c(0, main_y_max * 1.18), clip = "off") +
    ggplot2::labs(x = "Locked-anchor support", y = "Number of retained genes") +
    env$theme_imrs_publication(base_size = 12, legend_position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(lineheight = 0.95, margin = ggplot2::margin(t = 6)),
      plot.margin = ggplot2::margin(10, 14, 8, 10)
    )

  mini_y_max <- if (nrow(missing_anchor_tbl) > 0L) max(missing_anchor_tbl$n_retained_genes, na.rm = TRUE) else 1
  mini_plot <- ggplot2::ggplot(missing_anchor_tbl, ggplot2::aes(x = missing_anchor, y = n_retained_genes, fill = missing_anchor)) +
    ggplot2::geom_col(color = "#374151", linewidth = 0.35, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(n_retained_genes > 0, n_retained_genes, "")),
                       vjust = -0.25, size = 3.1, color = "#111827") +
    ggplot2::scale_fill_manual(values = env$anchor_palette[env$LOCKED_DATASETS_MOUSE], drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18)), breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::coord_cartesian(ylim = c(0, max(1, mini_y_max) * 1.22), clip = "off") +
    ggplot2::labs(x = "Missing anchor", y = "Retained genes") +
    env$theme_imrs_publication(base_size = 9.2, legend_position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 8.6),
      axis.title = ggplot2::element_text(size = 9.7),
      axis.title.x = ggplot2::element_text(size = 9.7, margin = ggplot2::margin(t = 6)),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  main_paths <- save_plot_all_formats_v5(main_plot, file.path(out_dir, "Figure2B_v5"), 5.8, 5.2, dpi)
  mini_paths <- save_plot_all_formats_v5(mini_plot, file.path(out_dir, "Figure2C_v5"), 3.6, 5.2, dpi)
  list(main = main_paths, mini = mini_paths)
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
        "1\u201324 h",
        old_env$format_time_label(time_h_num)
      ),
      display_label = dplyr::case_when(
        is_gse264344_acute ~ paste(
          "GSE264344 | adenoviral vector |",
          old_env$display_text(tissue),
          "| 1\u201324 h"
        ),
        is_gse264344_late ~ paste(
          "GSE264344 | adenoviral vector |",
          old_env$display_text(tissue),
          "|",
          old_env$format_time_label(time_h_num)
        ),
        TRUE ~ old_env$dataset_context_label(
          dataset_id, tissue, time_h, delivery_platform_clean, compact = TRUE
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
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 3)
      ),
      size = ggplot2::guide_legend(
        title = "Passing split contrasts",
        order = 2,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    old_env$theme_imrs_publication(base_size = 9.75) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 9.5),
      axis.title.x = ggplot2::element_text(size = 10.5),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.margin = ggplot2::margin(t = 4, r = 8, b = 4, l = 8),
      legend.box.margin = ggplot2::margin(t = 4, r = 4, b = 8, l = 4),
      plot.margin = ggplot2::margin(8, 18, 24, 12)
    )

  paths <- save_plot_all_formats_v5(
    strip_publication_text_v5(p),
    file.path(v5_root, "Figure1B_dataset_tissue_response_landscape_v5_corrected"),
    7.8, 6.8, dpi
  )

  log_msg_v5(
    "Corrected Figure 1B display table contributes ",
    nrow(dplyr::filter(plot_tbl, dataset_id == "GSE264344")),
    " GSE264344 rows."
  )
  list(paths = paths, plot_table = plot_tbl, plot_table_path = plot_table_path)
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
    height = 8.0,
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
  expected_panel_files <- unlist(lapply(
    specs$output_stem,
    function(stem) file.path(out_dir, paste0(stem, c(".png", ".pdf", ".svg")))
  ), use.names = FALSE)
  stale_panel_files <- expected_panel_files[file.exists(expected_panel_files)]
  if (length(stale_panel_files) > 0L && !all(file.remove(stale_panel_files))) {
    stop("Could not clear stale intermediate panel output(s) before rendering.", call. = FALSE)
  }
  lookup <- split(specs, specs$source_old_panel)
  rendered <- list()

  record <- function(spec, paths, width, height, generator) {
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
    record(spec, paths, spec$width, spec$height, source_code_section_or_function)
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
  record(spec_main, split_paths$main, spec_main$width, spec_main$height, "split_Figure2B_main")
  spec_mini <- specs[specs$source_function == "split_Figure2B_missing_anchor", , drop = FALSE]
  record(spec_mini, split_paths$mini, spec_mini$width, spec_mini$height, "split_Figure2B_missing_anchor")

  run_specs <- specs[!grepl("^split_Figure2B", specs$source_function), , drop = FALSE]
  for (i in seq_len(nrow(run_specs))) {
    fn_name <- run_specs$source_function[i]
    if (identical(fn_name, "render_Figure1A_merged_workflow_v5")) {
      log_msg_v5("Rendering v5 panel ", run_specs$panel_id[i],
                 " from merged workflow source ", run_specs$source_old_panel[i])
      paths <- render_figure1a_merged_workflow_v5(project_root, v5_root, dpi, workflow_script)
      record(run_specs[i, , drop = FALSE], paths,
             run_specs$width[i], run_specs$height[i], fn_name)
      next
    }
    if (identical(fn_name, "make_Figure1C_display_aggregated_v5")) {
      log_msg_v5("Rendering v5 panel ", run_specs$panel_id[i],
                 " from display-aggregated ", run_specs$source_old_panel[i])
      corrected <- render_corrected_figure1b_landscape_v5(old_env, v5_root, dpi)
      record(run_specs[i, , drop = FALSE], corrected$paths,
             run_specs$width[i], run_specs$height[i], fn_name)
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
  panel_manifest
}

assemble_v5_figures <- function(v5_root, panel_manifest, dpi = 400) {
  stop_if_missing_packages_v5(c("png", "grid", "svglite"))
  figure_plan <- make_v5_figure_plan()
  rows <- list()

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
        draw_png_fit_v5(hit$output_png[1], x, y, right - x, top - y)
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
      role = ifelse(grepl("^FigureS_", row$figure_id), "supplement", "main"),
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

validate_v5_outputs <- function(v5_root, manifest, baseline_v2, baseline_v3, baseline_v4) {
  expected_final <- c(
    "Figure1_main_v5", "Figure2_main_v5", "Figure3_main_v5",
    "Figure4_main_v5", "Figure5_main_v5",
    "FigureS_validation_faceted_summary_v5",
    "FigureS_weak_context_interpretation_categories_v5",
    "FigureS_comparator_benchmarking_v5"
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
