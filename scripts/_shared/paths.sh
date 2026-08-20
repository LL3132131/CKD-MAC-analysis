#!/usr/bin/env bash
# =============================================================================
# Title: paths.sh
# Purpose: Load and validate the shared release path configuration for shell.
# Manuscript Figure/Table: Infrastructure for Figures 1-6 and supplementary items
# Release role: SHARED INFRASTRUCTURE
# Required inputs: config/paths.env
# Generated outputs: Exported path variables; no analysis output
# Upstream dependencies: config/paths.env.example
# Downstream consumers: All public shell scripts
# Configuration keys: All keys in config/paths.env
# Expected environment: Bash
# Example run command: source scripts/_shared/paths.sh
# Reproducibility notes: This helper performs path resolution only.
# =============================================================================

load_project_paths() {
  local config_file="${PROJECT_PATHS_FILE:-config/paths.env}"
  if [[ ! -f "${config_file}" ]]; then
    echo "Path configuration not found: ${config_file}. Copy config/paths.env.example to config/paths.env and edit it." >&2
    return 2
  fi
  set -a
  # shellcheck disable=SC1090
  source "${config_file}"
  set +a

  local required=(
    PROJECT_ROOT DATA_RAW_DIR DATA_PROCESSED_DIR DATA_EXTERNAL_DIR RESULTS_DIR
    REFERENCE_DIR SCENIC_DB_DIR SCENIC_TF_LIST
    SCENIC_ALL_CELL_LOOM SCENIC_ALL_CELL_PREFIX
    SCENIC_SMC_LOOM SCENIC_SMC_PREFIX SCENIC_SMC_AUC_LOOM
    CREB3L2_ACTIVITY_CSV
    CUTTAG_FASTQ_DIR GWAS_DIR
    LD_REFERENCE_DIR MONOCLE2_PRIVATE_RDS MONOCLE2_PUBLIC_RDS
  )
  local key
  for key in "${required[@]}"; do
    if [[ -z "${!key:-}" ]]; then
      echo "Missing required configuration key: ${key}" >&2
      return 2
    fi
  done
}

require_input() {
  local input_path="$1"
  local label="${2:-input}"
  if [[ ! -e "${input_path}" ]]; then
    echo "Required ${label} does not exist: ${input_path}" >&2
    return 2
  fi
}

require_directory() {
  local directory_path="$1"
  local label="${2:-directory}"
  if [[ ! -d "${directory_path}" ]]; then
    echo "Required ${label} does not exist: ${directory_path}" >&2
    return 2
  fi
}
