# =============================================================================
# Title: paths.R
# Purpose: Load and validate the shared release path configuration for R.
# Manuscript Figure/Table: Infrastructure for Figures 1-6 and supplementary items
# Release role: SHARED INFRASTRUCTURE
# Required inputs: config/paths.env
# Generated outputs: Named path configuration in memory; no analysis output
# Upstream dependencies: config/paths.env.example
# Downstream consumers: All public R scripts
# Configuration keys: All keys in config/paths.env
# Expected environment: Base R
# Example run command: source("scripts/_shared/paths.R")
# Reproducibility notes: This helper performs path resolution only.
# =============================================================================

.expand_path_vars <- function(value, values) {
  for (i in seq_len(20L)) {
    hits <- regmatches(value, gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", value))[[1]]
    if (identical(hits, character(0)) || identical(hits, "-1")) break
    old <- value
    for (hit in unique(hits)) {
      key <- substring(hit, 3L, nchar(hit) - 1L)
      replacement <- values[[key]]
      if (is.null(replacement)) replacement <- Sys.getenv(key, unset = "")
      if (!nzchar(replacement)) {
        stop(sprintf("Undefined configuration key referenced by value: %s", key), call. = FALSE)
      }
      value <- sub(hit, replacement, value, fixed = TRUE)
    }
    if (identical(old, value)) break
  }
  value
}

load_project_paths <- function(config_file = Sys.getenv("PROJECT_PATHS_FILE", unset = "config/paths.env")) {
  if (!file.exists(config_file)) {
    stop(
      sprintf(
        "Path configuration not found: %s. Copy config/paths.env.example to config/paths.env and edit it.",
        config_file
      ),
      call. = FALSE
    )
  }
  lines <- readLines(config_file, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  values <- list()
  for (line in lines) {
    pos <- regexpr("=", line, fixed = TRUE)[1]
    if (pos < 2L) stop(sprintf("Invalid config line: %s", line), call. = FALSE)
    key <- trimws(substr(line, 1L, pos - 1L))
    value <- trimws(substr(line, pos + 1L, nchar(line)))
    value <- sub("^['\"]", "", value)
    value <- sub("['\"]$", "", value)
    values[[key]] <- value
  }
  for (key in names(values)) values[[key]] <- .expand_path_vars(values[[key]], values)

  required <- c(
    "PROJECT_ROOT", "DATA_RAW_DIR", "DATA_PROCESSED_DIR", "DATA_EXTERNAL_DIR",
    "RESULTS_DIR", "REFERENCE_DIR", "SCENIC_DB_DIR", "SCENIC_TF_LIST",
    "SCENIC_ALL_CELL_LOOM", "SCENIC_ALL_CELL_PREFIX",
    "SCENIC_SMC_LOOM", "SCENIC_SMC_PREFIX", "SCENIC_SMC_AUC_LOOM",
    "CREB3L2_ACTIVITY_CSV",
    "CUTTAG_FASTQ_DIR", "GWAS_DIR", "LD_REFERENCE_DIR",
    "MONOCLE2_PRIVATE_RDS", "MONOCLE2_PUBLIC_RDS"
  )
  missing <- required[!vapply(required, function(x) !is.null(values[[x]]) && nzchar(values[[x]]), logical(1))]
  if (length(missing)) stop(sprintf("Missing required configuration keys: %s", paste(missing, collapse = ", ")), call. = FALSE)

  config_dir <- dirname(normalizePath(config_file, winslash = "/", mustWork = TRUE))
  repo_root <- dirname(config_dir)
  if (!grepl("^([A-Za-z]:[/\\\\]|/)", values$PROJECT_ROOT)) {
    values$PROJECT_ROOT <- file.path(repo_root, values$PROJECT_ROOT)
  }
  values$PROJECT_ROOT <- normalizePath(values$PROJECT_ROOT, winslash = "/", mustWork = FALSE)
  for (key in names(values)) {
    if (grepl("(_DIR|_RDS|_LIST|_ROOT|_LOOM|_PREFIX|_CSV)$", key) && !grepl("^([A-Za-z]:[/\\\\]|/)", values[[key]])) {
      values[[key]] <- file.path(values$PROJECT_ROOT, values[[key]])
    }
  }
  values
}

configured_path <- function(paths, key, ..., trailing = FALSE) {
  base <- paths[[key]]
  if (is.null(base) || !nzchar(base)) stop(sprintf("Configuration key is undefined: %s", key), call. = FALSE)
  value <- file.path(base, ...)
  if (trailing) paste0(value, .Platform$file.sep) else value
}

require_input <- function(path, label = "input") {
  if (!file.exists(path)) stop(sprintf("Required %s does not exist: %s", label, path), call. = FALSE)
  path
}

require_directory <- function(path, label = "directory") {
  if (!dir.exists(path)) stop(sprintf("Required %s does not exist: %s", label, path), call. = FALSE)
  path
}

release_setwd <- function(path) {
  setwd(require_directory(path, "working/output directory"))
  invisible(path)
}
