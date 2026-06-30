# Shared path/config helpers for repository R scripts.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


locate_r_script_dir <- function(default = getwd()) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }

  script_file <- tryCatch(sys.frame(1)$ofile, error = function(...) NULL)
  if (!is.null(script_file)) {
    return(dirname(normalizePath(script_file, mustWork = FALSE)))
  }

  normalizePath(default, mustWork = FALSE)
}


find_project_root <- function(start_dir = NULL) {
  current <- normalizePath(start_dir %||% locate_r_script_dir(), mustWork = FALSE)
  for (i in seq_len(12)) {
    if (file.exists(file.path(current, "paths.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  stop("Could not find project root containing paths.yaml", call. = FALSE)
}


load_project_config <- function(project_root = NULL) {
  project_root <- project_root %||% find_project_root()
  config_file <- file.path(project_root, "paths.yaml")
  if (!file.exists(config_file)) {
    stop("Could not find root paths.yaml at ", config_file, call. = FALSE)
  }
  config <- yaml::read_yaml(config_file)
  attr(config, "config_file") <- config_file
  config
}


resolve_config_path <- function(project_root, config_value, default = NULL) {
  value <- config_value %||% default
  if (is.null(value)) {
    return(NULL)
  }
  if (is.character(value) && nzchar(value[[1]])) {
    return(file.path(project_root, value[[1]]))
  }
  if (is.list(value) && !is.null(value$root) && is.character(value$root) && nzchar(value$root[[1]])) {
    return(file.path(project_root, value$root[[1]]))
  }
  stop("Config value does not resolve to a path-like entry.", call. = FALSE)
}


prefer_active_conda_r_library <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
  if (!nzchar(conda_prefix)) {
    return(invisible(NULL))
  }

  conda_lib <- normalizePath(file.path(conda_prefix, "lib", "R", "library"), mustWork = FALSE)
  if (dir.exists(conda_lib)) {
    .libPaths(c(conda_lib, .libPaths()))
  }

  invisible(NULL)
}
