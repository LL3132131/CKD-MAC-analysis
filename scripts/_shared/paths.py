# =============================================================================
# Title: paths.py
# Purpose: Load and validate the shared release path configuration for Python.
# Manuscript Figure/Table: Infrastructure for Figures 1-6 and supplementary items
# Release role: SHARED INFRASTRUCTURE
# Required inputs: config/paths.env
# Generated outputs: Path configuration dictionary; no analysis output
# Upstream dependencies: config/paths.env.example
# Downstream consumers: All public Python scripts
# Configuration keys: All keys in config/paths.env
# Expected environment: Python 3 standard library
# Example run command: imported by release Python scripts
# Reproducibility notes: This helper performs path resolution only.
# =============================================================================
"""Shared path configuration for release Python scripts."""

from __future__ import annotations

import os
import re
from pathlib import Path

_VAR = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")
_REQUIRED = (
    "PROJECT_ROOT",
    "DATA_RAW_DIR",
    "DATA_PROCESSED_DIR",
    "DATA_EXTERNAL_DIR",
    "RESULTS_DIR",
    "REFERENCE_DIR",
    "SCENIC_DB_DIR",
    "SCENIC_TF_LIST",
    "SCENIC_ALL_CELL_LOOM",
    "SCENIC_ALL_CELL_PREFIX",
    "SCENIC_SMC_LOOM",
    "SCENIC_SMC_PREFIX",
    "SCENIC_SMC_AUC_LOOM",
    "CREB3L2_ACTIVITY_CSV",
    "CUTTAG_FASTQ_DIR",
    "GWAS_DIR",
    "LD_REFERENCE_DIR",
    "MONOCLE2_PRIVATE_RDS",
    "MONOCLE2_PUBLIC_RDS",
)


def load_project_paths(config_file: str | os.PathLike[str] | None = None) -> dict[str, str]:
    config_path = Path(config_file or os.environ.get("PROJECT_PATHS_FILE", "config/paths.env"))
    if not config_path.is_file():
        raise FileNotFoundError(
            f"Path configuration not found: {config_path}. "
            "Copy config/paths.env.example to config/paths.env and edit it."
        )
    values: dict[str, str] = {}
    for raw in config_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Invalid config line: {raw}")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip("'\"")

    def expand(value: str) -> str:
        for _ in range(20):
            updated = _VAR.sub(lambda m: values.get(m.group(1), os.environ.get(m.group(1), "")), value)
            if updated == value:
                break
            value = updated
        unresolved = _VAR.findall(value)
        if unresolved:
            raise KeyError(f"Undefined configuration keys: {', '.join(unresolved)}")
        return value

    values = {key: expand(value) for key, value in values.items()}
    missing = [key for key in _REQUIRED if not values.get(key)]
    if missing:
        raise KeyError(f"Missing required configuration keys: {', '.join(missing)}")

    repo_root = config_path.resolve().parent.parent
    project_root = Path(values["PROJECT_ROOT"])
    if not project_root.is_absolute():
        project_root = repo_root / project_root
    values["PROJECT_ROOT"] = str(project_root.resolve())
    for key, value in tuple(values.items()):
        if key.endswith(("_DIR", "_RDS", "_LIST", "_ROOT", "_LOOM", "_PREFIX", "_CSV")):
            path = Path(value)
            if not path.is_absolute():
                values[key] = str(project_root / path)
    return values


def configured_path(paths: dict[str, str], key: str, *parts: str) -> Path:
    if not paths.get(key):
        raise KeyError(f"Configuration key is undefined: {key}")
    return Path(paths[key]).joinpath(*parts)


def require_input(path: str | os.PathLike[str], label: str = "input") -> Path:
    value = Path(path)
    if not value.exists():
        raise FileNotFoundError(f"Required {label} does not exist: {value}")
    return value
