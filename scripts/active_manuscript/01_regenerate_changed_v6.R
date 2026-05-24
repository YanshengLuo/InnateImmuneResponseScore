#!/usr/bin/env Rscript

# Regenerate only v6 figures with plot-visible text changes.

options(stringsAsFactors = FALSE)

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
                      error = function(e) NA_character_)
if (is.na(this_file)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
}
active_config_helper <- file.path(dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE)),
                                  "_active_config.R")
if (!file.exists(active_config_helper)) {
  active_config_helper <- file.path(getwd(), "scripts", "active_manuscript", "_active_config.R")
}
source(active_config_helper)

config <- imrs_load_active_config(dirname(active_config_helper))
project_root <- imrs_project_root(config)
v6_root <- imrs_config_field_path(config, "manuscript_output_dir")
helper <- imrs_config_field_path(config, "figure_helper_script", "scripts/utilities/v5_helpers.R")

source(helper)
stop_if_missing_packages_v5()

dpi <- 400
old_script <- imrs_config_field_path(config, "legacy_figure_generator_script")
if (!file.exists(old_script)) {
  stop("Legacy figure generator is required for this targeted regeneration step but was not found: ",
       old_script, call. = FALSE)
}
old_env <- load_old_generator_env_v5(old_script, project_root, v6_root)
specs <- make_v5_panel_plan()
out_dir <- file.path(v6_root, "intermediate_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_base_dir <- function(panel_id) {
  if (panel_id %in% c("Figure1A", "Figure1B")) v6_root else out_dir
}

paths_for_spec <- function(spec) {
  base <- file.path(panel_base_dir(spec$panel_id), spec$output_stem)
  c(png = paste0(base, ".png"), pdf = paste0(base, ".pdf"), svg = paste0(base, ".svg"))
}

panel_manifest_from_specs <- function() {
  rows <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, , drop = FALSE]
    paths <- paths_for_spec(spec)
    rows[[i]] <- data.frame(
      figure_id = spec$figure_id,
      panel_id = spec$panel_id,
      panel_letter = spec$panel_letter,
      source_old_panel = spec$source_old_panel,
      source_function = spec$source_function,
      output_png = norm_path_v5(paths[["png"]], must_work = TRUE),
      output_pdf = norm_path_v5(paths[["pdf"]], must_work = TRUE),
      output_svg = norm_path_v5(paths[["svg"]], must_work = TRUE),
      width = spec$width,
      height = spec$height,
      dpi = dpi,
      notes = spec$notes,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

panel_manifest <- panel_manifest_from_specs()
regenerated_panels <- character()

record_panel <- function(spec, paths, generator) {
  idx <- match(spec$panel_id, panel_manifest$panel_id)
  panel_manifest[idx, "source_function"] <<- generator
  panel_manifest[idx, "output_png"] <<- norm_path_v5(paths[["png"]], must_work = TRUE)
  panel_manifest[idx, "output_pdf"] <<- norm_path_v5(paths[["pdf"]], must_work = TRUE)
  panel_manifest[idx, "output_svg"] <<- norm_path_v5(paths[["svg"]], must_work = TRUE)
  regenerated_panels <<- unique(c(regenerated_panels, spec$panel_id))
}

target_specs <- specs[specs$panel_id %in% c("Figure3B", "Figure5A"), , drop = FALSE]
lookup <- split(target_specs, target_specs$source_old_panel)

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
  record_panel(spec, paths, source_code_section_or_function)
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
  paths <- save_grid_all_formats_v5(out_base, spec$width, spec$height, dpi, draw_fun)
  record_panel(spec, paths, source_code_section_or_function)
  invisible(paths)
}

spec_1a <- specs[specs$panel_id == "Figure1A", , drop = FALSE]
paths_1a <- render_figure1a_merged_workflow_v5(project_root, v6_root, dpi)
record_panel(spec_1a, paths_1a, "render_Figure1A_merged_workflow_v5")

spec_1b <- specs[specs$panel_id == "Figure1B", , drop = FALSE]
paths_1b <- render_corrected_figure1b_landscape_v5(old_env, v6_root, dpi)$paths
record_panel(spec_1b, paths_1b, "make_Figure1C_display_aggregated_v5")

old_env$make_FigureSB_simplified()
old_env$make_Figure4A()

write_tsv_v5(panel_manifest, file.path(v6_root, "tables", "v5_panel_manifest.tsv"))

figure_manifest_from_plan <- function() {
  figure_plan <- make_v5_figure_plan()
  rows <- lapply(figure_plan, function(fig) {
    data.frame(
      figure_id = fig$figure_id,
      output_stem = fig$output_stem,
      role = fig$role,
      output_png = norm_path_v5(file.path(v6_root, paste0(fig$output_stem, ".png")), must_work = TRUE),
      output_pdf = norm_path_v5(file.path(v6_root, paste0(fig$output_stem, ".pdf")), must_work = TRUE),
      output_svg = norm_path_v5(file.path(v6_root, paste0(fig$output_stem, ".svg")), must_work = TRUE),
      width = fig$width,
      height = fig$height,
      dpi = dpi,
      notes = "Combined v6 figure contains clean panel graphics only; no large internal plot titles.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

assemble_selected_figures <- function(selected_figures) {
  figure_plan <- make_v5_figure_plan()
  figure_plan <- figure_plan[vapply(figure_plan, function(fig) fig$figure_id %in% selected_figures, logical(1))]
  rows <- list()

  for (fig in figure_plan) {
    width <- fig$width
    height <- fig$height
    png_path <- file.path(v6_root, paste0(fig$output_stem, ".png"))
    pdf_path <- file.path(v6_root, paste0(fig$output_stem, ".pdf"))
    svg_path <- file.path(v6_root, paste0(fig$output_stem, ".svg"))

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
          stop("Missing v6 panel PNG for ", panel_id, call. = FALSE)
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
    save_grid_device_v5(svg_path, width, height, dpi, draw_fun, "svg")

    rows[[length(rows) + 1L]] <- data.frame(
      figure_id = fig$figure_id,
      output_stem = fig$output_stem,
      role = fig$role,
      output_png = norm_path_v5(png_path, must_work = TRUE),
      output_pdf = norm_path_v5(pdf_path, must_work = TRUE),
      output_svg = norm_path_v5(svg_path, must_work = TRUE),
      width = width,
      height = height,
      dpi = dpi,
      notes = "Combined v6 figure contains clean panel graphics only; no large internal plot titles.",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

figure_manifest <- figure_manifest_from_plan()
selected_rows <- assemble_selected_figures(c("Figure1", "Figure3", "Figure5"))
figure_manifest[match(selected_rows$figure_id, figure_manifest$figure_id), ] <- selected_rows
write_tsv_v5(figure_manifest, file.path(v6_root, "tables", "v5_figure_manifest_wide.tsv"))

manifest <- manifest_long_v5(panel_manifest, figure_manifest, v6_root)
file.copy(file.path(v6_root, "figure_v5_manifest.tsv"),
          file.path(v6_root, "figure_v6_manifest.tsv"),
          overwrite = TRUE)

log_lines <- c(
  paste0("start_time\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("source_folder\t", imrs_config_field_path(config, "source_v5_dir")),
  paste0("output_folder\t", v6_root),
  paste0("regenerated_panels\t", paste(regenerated_panels, collapse = ", ")),
  paste0("regenerated_figures\t", paste(c("Figure1", "Figure3", "Figure5"), collapse = ", ")),
  paste0("manifest_records\t", nrow(manifest))
)
writeLines(log_lines, file.path(v6_root, "v6_regeneration_targeted_run.log"), useBytes = TRUE)

cat("Regenerated v6 panels: ", paste(regenerated_panels, collapse = ", "), "\n", sep = "")
cat("Regenerated v6 figures: Figure1, Figure3, Figure5\n")
