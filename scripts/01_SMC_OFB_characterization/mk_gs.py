# =============================================================================
# Title: mk_gs.py
# Purpose: scDRS gene-set munging
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: MAGMA gene results and gene locations
# Generated outputs: MAGMA Z-score tables with symbols
# Upstream dependencies: MAGMA gene results and gene locations
# Downstream consumers: scDRS gene-set munging
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/01_SMC_OFB_characterization/mk_gs.py
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_shared = _ReleasePath(__file__).resolve().parents[1] / "_shared"
_release_sys.path.insert(0, str(_release_shared))
from paths import configured_path, load_project_paths
PATHS = load_project_paths()
import pandas as pd

step2_path = configured_path(
    PATHS,
    "SMC_RESULTS_DIR",
    "magma",
    "step2.genes.out",
)
if not step2_path.is_file() or step2_path.stat().st_size == 0:
    raise FileNotFoundError(f"MAGMA gene result is missing or empty: {step2_path}")

magma_input = pd.read_csv(str(step2_path), sep=r"\s+")
required_magma_columns = {"GENE", "ZSTAT"}
missing_magma_columns = sorted(required_magma_columns - set(magma_input.columns))
if missing_magma_columns:
    raise ValueError(
        "MAGMA gene result is missing required columns: "
        + ", ".join(missing_magma_columns)
    )

magma_df = magma_input.loc[:, ["GENE", "ZSTAT"]].rename(
    columns={"ZSTAT": "zscore"}
)
magma_df = magma_df.dropna(subset=["zscore"])
if not pd.api.types.is_numeric_dtype(magma_df["zscore"]):
    raise TypeError("MAGMA ZSTAT values must be numeric.")

magma_df.to_csv(str(configured_path(PATHS, "SMC_RESULTS_DIR", "magma", "magma_zscore.tsv")), sep="\t", index=False)
######################注释基因符号##############################
# 读取NCBI37.3.gene.loc
gene_loc_df = pd.read_csv(
    str(configured_path(PATHS, "GWAS_DIR", "NCBI38", "NCBI38.gene.loc")),
    sep=r"\s+",  # 通用空白符分隔
    header=None,  # 无列名
    usecols=[0, 5],  # 提取第一列（Entrez ID）和最后一列（基因符号）
    names=["entrez_id", "symbol"]  # 重命名列
)

# 建立Entrez ID到基因符号的字典（去除重复或无效项）
entrez_to_symbol = dict(zip(gene_loc_df["entrez_id"], gene_loc_df["symbol"]))

# 读取原始magma_zscore.tsv（包含Entrez ID）
magma_df = pd.read_csv(
    str(configured_path(PATHS, "SMC_RESULTS_DIR", "magma", "magma_zscore.tsv")),
    sep="\t"
)
# 将Entrez ID转换为基因符号	
magma_df["GENE"] = magma_df["GENE"].map(entrez_to_symbol)

magma_df = magma_df.dropna(subset=["GENE"])
if magma_df.empty:
    raise ValueError("No MAGMA genes could be mapped to gene symbols.")
magma_df["GENE"] = magma_df["GENE"].astype(str)
if magma_df["GENE"].isna().any() or magma_df["GENE"].str.strip().eq("").any():
    raise ValueError("Mapped GENE values must not be missing or empty.")
if not pd.api.types.is_numeric_dtype(magma_df["zscore"]):
    raise TypeError("Mapped MAGMA zscore values must be numeric.")
if list(magma_df.columns) != ["GENE", "zscore"]:
    raise ValueError("The final MAGMA table must contain exactly GENE and zscore.")
magma_df.to_csv(
    str(configured_path(PATHS, "SMC_RESULTS_DIR", "magma", "magma_zscore_symbol.tsv")),
    sep="\t",
    index=False
)
