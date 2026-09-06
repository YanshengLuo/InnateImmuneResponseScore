# Shared config helpers for the repository active-manuscript layer.

imrs_repo_root <- function(start = getwd()) {
  env_root <- Sys.getenv("IMRS_REPOSITORY_ROOT", unset = "")
  starts <- unique(c(env_root, start))
  starts <- starts[!is.na(starts) & nzchar(starts)]

  for (item in starts) {
    current <- normalizePath(item, winslash = "/", mustWork = FALSE)
    if (file.exists(current) && !dir.exists(current)) current <- dirname(current)
    repeat {
      if (file.exists(file.path(current, "config", "config_template.yml")) &&
          file.exists(file.path(current, "scripts", "active_manuscript", "lib", "active_config.R"))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  stop(
    "Could not identify the IMRS repository root from: ", paste(starts, collapse = "; "),
    ". Keep the extracted repository structure intact.",
    call. = FALSE
  )
}

imrs_read_simple_yaml <- function(file) {
  if (!file.exists(file)) {
    stop("Config file does not exist: ", file, call. = FALSE)
  }
  lines <- readLines(file, warn = FALSE)
  lines <- lines[grepl(":", lines) & !grepl("^\\s*#", lines)]
  out <- list()
  for (line in lines) {
    key <- trimws(sub(":.*$", "", line))
    value <- trimws(sub("^[^:]+:", "", line))
    value <- trimws(gsub("^['\"]|['\"]$", "", value))
    if (tolower(value) %in% c("true", "false", "yes", "no")) {
      value <- tolower(value) %in% c("true", "yes")
    }
    out[[key]] <- value
  }
  out
}

imrs_load_active_config <- function(start = getwd()) {
  repo_root <- imrs_repo_root(start)
  config_override <- Sys.getenv("IMRS_ACTIVE_CONFIG", unset = "")
  if (nzchar(config_override)) {
    config_file <- normalizePath(config_override, winslash = "/", mustWork = FALSE)
  } else {
    config_file <- file.path(repo_root, "config", "config.yml")
    if (!file.exists(config_file)) {
      config_file <- file.path(repo_root, "config", "config_template.yml")
    }
  }
  config <- imrs_read_simple_yaml(config_file)
  config$repo_root <- repo_root
  config$config_file <- normalizePath(config_file, winslash = "/", mustWork = FALSE)
  if (is.null(config$project_root) || !nzchar(config$project_root)) {
    config$project_root <- "."
  }
  config
}

imrs_config_path <- function(config, value, default = NULL, must_work = FALSE) {
  if (is.null(value) || is.na(value) || !nzchar(value)) value <- default
  if (is.null(value) || is.na(value) || !nzchar(value)) return(NA_character_)
  base <- config$repo_root
  out <- if (grepl("^([A-Za-z]:|/)", value)) value else file.path(base, value)
  normalizePath(out, winslash = "/", mustWork = must_work)
}

imrs_config_field_path <- function(config, field, default = NULL, must_work = FALSE) {
  imrs_config_path(config, config[[field]], default = default, must_work = must_work)
}

imrs_project_root <- function(config) {
  imrs_config_path(config, config$project_root, default = ".", must_work = FALSE)
}

imrs_bool <- function(x) isTRUE(x) || identical(tolower(as.character(x)), "true")

imrs_first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}
