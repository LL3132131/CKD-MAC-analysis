# =============================================================================
# Title: Prepare processed RDS files for the public workflow
# Purpose: Map public RDS filenames to the local filenames expected by code
# Manuscript Figure/Table: Figures 1-6 (file preparation only)
# Release role: PUBLIC DATA PREPARATION HELPER
# Required inputs: One command-line directory containing the four public RDS files
# Generated outputs: Byte-identical copies under data/processed/ and a mapping log
# Upstream dependencies: Public processed RDS files from the finalized release location
# Downstream consumers: Existing public analysis and plotting scripts
# Configuration keys: None; the input directory is the first command-line argument
# Expected environment: Base R plus an operating-system SHA-256 utility
# Example run command:
#   Rscript scripts/00_prepare_processed_data.R path/to/downloaded_files
# Reproducibility notes: This script never reads or modifies RDS contents. It
# performs existence, size, byte-copy, and SHA-256 identity checks only.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop(
    "Usage: Rscript scripts/00_prepare_processed_data.R ",
    "<directory-containing-public-RDS-files>"
  )
}

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Unable to determine this script's location from the Rscript command.")
}

script_path <- normalizePath(
  sub("^--file=", "", script_arg),
  winslash = "/",
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = TRUE
)
source_dir <- normalizePath(
  args[[1]],
  winslash = "/",
  mustWork = TRUE
)
if (!isTRUE(file.info(source_dir)$isdir)) {
  stop("The command-line input must be a directory: ", source_dir)
}

sha256_file <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)

  commands <- list(
    list(command = Sys.which("sha256sum"), args = c(shQuote(path))),
    list(command = Sys.which("shasum"), args = c("-a", "256", shQuote(path))),
    list(command = Sys.which("certutil"), args = c("-hashfile", shQuote(path), "SHA256"))
  )

  for (candidate in commands) {
    if (!nzchar(candidate$command)) {
      next
    }
    output <- suppressWarnings(
      system2(
        candidate$command,
        args = candidate$args,
        stdout = TRUE,
        stderr = TRUE
      )
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    if (!identical(as.integer(status), 0L)) {
      next
    }
    hashes <- unlist(
      regmatches(
        output,
        gregexpr("(?i)\\b[0-9a-f]{64}\\b", output, perl = TRUE)
      ),
      use.names = FALSE
    )
    hashes <- unique(tolower(hashes[nzchar(hashes)]))
    if (length(hashes) == 1L) {
      return(hashes)
    }
  }

  stop(
    "No supported SHA-256 utility was available. Install or expose one of ",
    "sha256sum, shasum, or certutil; no file was copied."
  )
}

mapping <- data.frame(
  public_filename = c(
    "8data_cca_3000_public.rds",
    "7data_umap_public.rds",
    "smc_subclusters_final_seurat_public.rds",
    "sub3_m2_public.rds"
  ),
  local_code_filename = c(
    "8data_cca_3000.rds",
    "7data_umap.rds",
    "smc_subclusters_final_seurat.rds",
    "sub3_m2_public.rds"
  ),
  stringsAsFactors = FALSE
)

mapping$source_path <- file.path(source_dir, mapping$public_filename)
target_dir <- file.path(project_root, "data", "processed")
mapping$target_path <- file.path(target_dir, mapping$local_code_filename)

missing_sources <- mapping$source_path[
  !file.exists(mapping$source_path) |
    is.na(file.info(mapping$source_path)$size) |
    file.info(mapping$source_path)$size <= 0
]
if (length(missing_sources) > 0L) {
  stop(
    "Required public RDS files are missing or empty: ",
    paste(basename(missing_sources), collapse = ", ")
  )
}

mapping$size_bytes <- as.numeric(file.info(mapping$source_path)$size)
mapping$sha256 <- vapply(mapping$source_path, sha256_file, character(1))

# Complete the conflict preflight before creating or copying anything.
for (index in seq_len(nrow(mapping))) {
  target <- mapping$target_path[[index]]
  if (file.exists(target)) {
    if (isTRUE(file.info(target)$isdir)) {
      stop("Expected a file but found a directory: ", target)
    }
    target_hash <- sha256_file(target)
    if (!identical(target_hash, mapping$sha256[[index]])) {
      stop(
        "Refusing to overwrite an existing file with a different SHA-256: ",
        target
      )
    }
  }
}

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(target_dir)) {
  stop("Failed to create output directory: ", target_dir)
}

for (index in seq_len(nrow(mapping))) {
  source <- mapping$source_path[[index]]
  target <- mapping$target_path[[index]]

  if (!file.exists(target)) {
    temporary_target <- paste0(target, ".partial.", Sys.getpid())
    if (file.exists(temporary_target)) {
      stop("Temporary copy target already exists: ", temporary_target)
    }
    copied <- file.copy(
      from = source,
      to = temporary_target,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )
    if (!isTRUE(copied)) {
      stop("Byte-copy failed for: ", basename(source))
    }
    temporary_hash <- sha256_file(temporary_target)
    if (!identical(temporary_hash, mapping$sha256[[index]])) {
      unlink(temporary_target)
      stop("SHA-256 mismatch after byte-copy: ", basename(source))
    }
    if (file.exists(target)) {
      unlink(temporary_target)
      stop("Target appeared during copying; refusing to overwrite: ", target)
    }
    if (!isTRUE(file.rename(temporary_target, target))) {
      unlink(temporary_target)
      stop("Failed to finalize copied file: ", target)
    }
  }

  final_hash <- sha256_file(target)
  final_size <- as.numeric(file.info(target)$size)
  if (!identical(final_hash, mapping$sha256[[index]]) ||
      !identical(final_size, mapping$size_bytes[[index]])) {
    stop("Final byte-level validation failed for: ", target)
  }
}

public_log <- data.frame(
  public_filename = mapping$public_filename,
  local_code_filename = mapping$local_code_filename,
  size_bytes = mapping$size_bytes,
  sha256 = mapping$sha256,
  status = "VERIFIED_BYTE_IDENTICAL",
  stringsAsFactors = FALSE
)

log_path <- file.path(target_dir, "public_filename_mapping.tsv")
temporary_log <- paste0(log_path, ".partial.", Sys.getpid())
if (file.exists(temporary_log)) {
  stop("Temporary log target already exists: ", temporary_log)
}
write.table(
  public_log,
  file = temporary_log,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  fileEncoding = "UTF-8"
)

if (file.exists(log_path)) {
  if (!identical(sha256_file(log_path), sha256_file(temporary_log))) {
    unlink(temporary_log)
    stop("Existing mapping log differs; refusing to overwrite: ", log_path)
  }
  unlink(temporary_log)
} else if (!isTRUE(file.rename(temporary_log, log_path))) {
  unlink(temporary_log)
  stop("Failed to finalize mapping log: ", log_path)
}

message(
  "Verified four byte-identical public RDS mappings. Log: ",
  normalizePath(log_path, winslash = "/", mustWork = TRUE)
)
