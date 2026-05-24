build_merged_imrs_workflow_v5 <- function(project_root = normalizePath(Sys.getenv("IMRS_PROJECT_ROOT", "."),
                                                                       winslash = "/", mustWork = FALSE),
                                          output_dir = Sys.getenv("IMRS_MANUSCRIPT_OUTPUT_DIR",
                                                                  file.path(project_root, "results_release_templates")),
                                          output_stem = "Figure1A_IMRS_merged_workflow_v5",
                                          width = 7.5,
                                          height = 8.0,
                                          dpi = 400) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to build the merged IMRS workflow.", call. = FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required to build the merged IMRS workflow.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  box_text_size <- 3.2
  column_text_size <- 2.95
  connector_text_size <- 2.8
  title_text_size <- 4.25
  subtitle_text_size <- 2.8
  footer_text_size <- 2.85

  node_border <- "#24364A"
  arrow_color <- "#46556A"
  text_color <- "#111827"

  nodes <- data.frame(
    id = c(
      "inputs", "split",
      "anchor_de", "anchor_reproducible", "anchor_filter", "frozen_coefficients",
      "target_norm", "target_z", "target_raw", "target_imrsz", "target_eval",
      "role_audit", "benchmarking", "interpretation"
    ),
    label = c(
      "Verified metadata and raw\nRNA-seq count matrices",
      "Delivery-versus-control\nsplit definitions",
      "Locked-anchor\ndelivery-versus-control\ndifferential expression",
      "Cross-anchor\nreproducibility filtering",
      "Heterogeneity and\ninformation-content filtering",
      "Frozen anchor-derived\ngene coefficients",
      "Target-dataset\ncount normalization",
      "Control-referenced\ngene z-scores",
      "Weighted sample-level\nIMRS score",
      "Control-standardized\nsample IMRSz",
      "Delivery-minus-control\nΔIMRSz, directionality,\nand secondary AUC",
      "Dataset-role curation and\nboundary-context audit",
      "Label permutation, comparator benchmarking,\nand coefficient sensitivity",
      "Manuscript\nanalysis groups"
    ),
    xmin = c(1.70, 1.70, 0.20, 0.20, 0.20, 0.20, 8.55, 8.55, 8.55, 8.55, 8.55, 1.55, 1.55, 1.55),
    xmax = c(12.30, 12.30, 5.35, 5.35, 5.35, 5.35, 13.80, 13.80, 13.80, 13.80, 13.80, 12.45, 12.45, 12.45),
    ymin = c(11.60, 10.55, 8.45, 7.30, 6.15, 5.00, 8.45, 7.30, 6.15, 5.00, 3.65, 2.80, 1.70, 0.65),
    ymax = c(12.25, 11.20, 9.35, 8.00, 6.85, 5.70, 9.35, 8.00, 6.85, 5.70, 4.65, 3.50, 2.45, 1.35),
    fill = c(
      "#EAF2F8", "#EAF2F8",
      "#E8F5E9", "#E8F5E9", "#E8F5E9", "#E8F5E9",
      "#FFF3E0", "#FFF3E0", "#FFF3E0", "#FFF3E0", "#F3E8FF",
      "#F8E8EE", "#F8E8EE", "#F8E8EE"
    ),
    stringsAsFactors = FALSE
  )
  nodes$x <- (nodes$xmin + nodes$xmax) / 2
  nodes$y <- (nodes$ymin + nodes$ymax) / 2

  arrows <- data.frame(
    x = c(
      7.00, 7.00,
      7.00, 7.00,
      2.78, 2.78, 2.78,
      11.18, 11.18, 11.18, 11.18,
      5.35,
      11.18, 7.00, 7.00
    ),
    y = c(
      11.60, 10.55,
      10.55, 10.55,
      8.45, 7.30, 6.15,
      8.45, 7.30, 6.15, 5.00,
      5.35,
      3.65, 2.80, 1.70
    ),
    xend = c(
      7.00, 7.00,
      2.78, 11.18,
      2.78, 2.78, 2.78,
      11.18, 11.18, 11.18, 11.18,
      8.55,
      7.00, 7.00, 7.00
    ),
    yend = c(
      11.20, 11.20,
      9.35, 9.35,
      8.00, 6.85, 5.70,
      8.00, 6.85, 5.70, 4.65,
      6.50,
      3.50, 2.45, 1.35
    ),
    stringsAsFactors = FALSE
  )

  arrow_spec <- grid::arrow(
    length = grid::unit(0.12, "inches"),
    type = "closed"
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.45,
      lineend = "round",
      color = arrow_color,
      arrow = arrow_spec
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      color = node_border,
      linewidth = 0.45
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label),
      size = box_text_size,
      lineheight = 0.90,
      color = text_color,
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 7.0, y = 14.05,
      label = "Frozen IMRS scoring and transfer-evaluation workflow",
      size = title_text_size,
      fontface = "bold",
      lineheight = 0.92,
      color = text_color,
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 7.0, y = 13.40,
      label = "Anchor-derived gene coefficients are frozen before scoring\nindependent delivery-versus-control contrasts",
      size = subtitle_text_size,
      lineheight = 0.96,
      color = "#374151",
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 1.80, y = 9.96,
      label = "Locked-anchor model\nconstruction",
      size = column_text_size,
      fontface = "bold",
      lineheight = 0.92,
      color = text_color,
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 12.05, y = 9.96,
      label = "Frozen-coefficient scoring\nand transfer evaluation",
      size = column_text_size,
      lineheight = 0.92,
      fontface = "bold",
      color = text_color,
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 7.0, y = 6.10,
      label = "frozen coefficients",
      size = connector_text_size,
      angle = 15,
      color = "#4B5563",
      family = "sans"
    ) +
    ggplot2::annotate(
      "text", x = 7.0, y = 0.22,
      label = "Frozen anchor-derived coefficients are not refit during validation or transfer evaluation.",
      size = footer_text_size,
      color = "#374151",
      family = "sans"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-0.15, 14.15),
      ylim = c(-0.05, 14.35),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::theme_void(base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(18, 20, 16, 20)
    )

  output_base <- file.path(output_dir, output_stem)
  paths <- c(
    png = paste0(output_base, ".png"),
    pdf = paste0(output_base, ".pdf"),
    svg = paste0(output_base, ".svg")
  )

  ggplot2::ggsave(paths[["png"]], p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paths[["pdf"]], p, width = width, height = height, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paths[["svg"]], p, width = width, height = height, device = svglite::svglite, bg = "white", limitsize = FALSE)
  } else {
    paths[["svg"]] <- NA_character_
  }

  paths
}
