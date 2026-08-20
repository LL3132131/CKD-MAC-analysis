# =============================================================================
# Title: plot.py
# Purpose: Figure 2
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: AnnData, gene set, score and group files
# Generated outputs: scDRS UMAP/group-statistic figures, including a focused
# OFB-SMC/Fibroblast/non-OFB-SMC UMAP
# Upstream dependencies: AnnData, gene set, score and group files
# Downstream consumers: Figure 2
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/01_SMC_OFB_characterization/plot.py
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_shared = _ReleasePath(__file__).resolve().parents[1] / "_shared"
_release_sys.path.insert(0, str(_release_shared))
from paths import configured_path, load_project_paths
PATHS = load_project_paths()
import scdrs
import scanpy as sc
sc.set_figure_params(dpi=300)
from anndata import AnnData
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings

warnings.filterwarnings("ignore")

outpath=str(configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs"))
os.chdir(outpath)

adata = sc.read_h5ad(str(configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs", "7data.h5ad")))
df_gs = pd.read_csv(str(configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs", "disease_geneset.gs")), sep="\t", index_col=0, encoding="utf-8")
traits = list(df_gs.index)
if not traits:
    raise ValueError("No traits were found in disease_geneset.gs.")
if pd.Index(traits).has_duplicates:
    raise ValueError("Trait names in disease_geneset.gs must be unique.")

if "Celltype" not in adata.obs.columns:
    raise ValueError("Required AnnData metadata column is missing: Celltype")
if adata.obs_names.has_duplicates:
    raise ValueError("AnnData cell barcodes must be unique.")

atlas_celltype_order = [
    "Neutrophil",
    "Fibroblast",
    "non-OFB-SMC",
    "Endothelial cell",
    "CD8 T cell",
    "OFB-SMC",
    "CD4 T cell",
    "NK Cell",
    "Pericyte",
    "Mixed myeloid cell",
    "Macrophage",
    "ProMacs",
    "B cell",
    "Eosinophil",
    "Mast cell",
]
prohibited_7data_celltypes = {"DC", "Dendritic cell"}
available_adata_celltypes = set(adata.obs["Celltype"].dropna().astype(str))
present_prohibited_celltypes = sorted(
    available_adata_celltypes & prohibited_7data_celltypes
)
if present_prohibited_celltypes:
    raise ValueError(
        "The scDRS AnnData uses an eight-sample or obsolete seven-sample "
        "Celltype interface: "
        + ", ".join(present_prohibited_celltypes)
    )
missing_adata_celltypes = [
    celltype
    for celltype in atlas_celltype_order
    if celltype not in available_adata_celltypes
]
unexpected_adata_celltypes = sorted(
    available_adata_celltypes - set(atlas_celltype_order)
)
if missing_adata_celltypes or unexpected_adata_celltypes:
    details = []
    if missing_adata_celltypes:
        details.append("missing: " + ", ".join(missing_adata_celltypes))
    if unexpected_adata_celltypes:
        details.append("unexpected: " + ", ".join(unexpected_adata_celltypes))
    raise ValueError(
        "AnnData Celltype values do not match the final seven-sample interface ("
        + "; ".join(details)
        + ")"
    )

dict_score = {
    trait: pd.read_csv(
        str(
            configured_path(
                PATHS,
                "SMC_RESULTS_DIR",
                "scdrs",
                f"{trait}.full_score.gz",
            )
        ),
        sep="\t",
        index_col=0,
        encoding="utf-8",
    )
    for trait in traits
}

for trait in dict_score:
    score_table = dict_score[trait]
    if "norm_score" not in score_table.columns:
        raise ValueError(f"Required scDRS score column is missing for {trait}: norm_score")
    if score_table.index.has_duplicates:
        raise ValueError(f"scDRS score table contains duplicate cell barcodes for {trait}.")
    missing_score_cells = adata.obs_names.difference(score_table.index)
    extra_score_cells = score_table.index.difference(adata.obs_names)
    if len(missing_score_cells) > 0 or len(extra_score_cells) > 0:
        raise ValueError(
            f"scDRS score cells do not exactly match AnnData cells for {trait} "
            f"(missing scores: {len(missing_score_cells)}; "
            f"extra scores: {len(extra_score_cells)})."
        )
    score_table = score_table.reindex(adata.obs_names)
    if score_table["norm_score"].isna().any():
        raise ValueError(f"scDRS norm_score contains missing values for {trait}.")
    dict_score[trait] = score_table
    adata.obs[trait] = score_table["norm_score"]

sc.set_figure_params(figsize=[2.5, 2.5], dpi=300)

sc.pl.umap(
    adata,
    color=dict_score.keys(),
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    s=20,
    show=False,
)
plt.savefig(str(configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs", "umap_plot.png")), dpi=300, bbox_inches='tight')
plt.close()
#plt.savefig(str(configured_path(PATHS, "RESULTS_DIR", "figures", "Figure2", "scDRS_umap_plot.pdf")), dpi=900, bbox_inches='tight')
##############
data_dir = str(configured_path(PATHS, "SMC_RESULTS_DIR", "scdrs"))
dict_df_stats = {
    trait: pd.read_csv(
        os.path.join(data_dir, f"{trait}.scdrs_group.Celltype"),
        sep="\t", 
        index_col=0,
        encoding="utf-8"
    )
    for trait in traits
}
# ------------ 2. 最终15类major-cell顺序与显式历史兼容映射 ------------
# Public inputs are expected to use atlas_celltype_order. This dictionary is
# retained only as an explicit, auditable compatibility path for historical
# scDRS group tables; it is never applied when the input already uses final names.
legacy_to_final = {
    "NE": "Neutrophil",
    "FB": "Fibroblast",
    "SMC": "non-OFB-SMC",
    "EC": "Endothelial cell",
    "CD8_T_Cell": "CD8 T cell",
    "CD4_help_T_Cell": "CD4 T cell",
    "NK": "NK Cell",
    "unKnown": "Pericyte",
    "Pericyte Cell": "Pericyte",
    "FB_SMC": "OFB-SMC",
    "Macrophage_NE": "Mixed myeloid cell",
    "Proliferation": "ProMacs",
    "TRMs": "ProMacs",
    "B_Cell": "B cell",
    "Eosinophils": "Eosinophil",
    "Mastocyte": "Mast cell",
}


def standardize_group_stats_index(df_stats, trait):
    required_stats_columns = {"assoc_mcz", "assoc_mcp"}
    missing_stats_columns = sorted(required_stats_columns - set(df_stats.columns))
    if missing_stats_columns:
        raise ValueError(
            f"Required scDRS group-statistic columns are missing for {trait}: "
            + ", ".join(missing_stats_columns)
        )

    standardized = df_stats.copy()
    original_index = standardized.index.astype(str)
    duplicate_original_labels = sorted(
        set(original_index[original_index.duplicated(keep=False)])
    )
    if duplicate_original_labels:
        raise ValueError(
            f"scDRS group table contains duplicate Celltype rows for {trait}: "
            + ", ".join(duplicate_original_labels)
        )
    present_prohibited_labels = sorted(
        set(original_index) & prohibited_7data_celltypes
    )
    if present_prohibited_labels:
        raise ValueError(
            f"scDRS group table uses an eight-sample or obsolete seven-sample "
            f"Celltype interface for {trait}: "
            + ", ".join(present_prohibited_labels)
        )
    pericyte_aliases = {"unKnown", "Pericyte Cell"}
    present_pericyte_aliases = sorted(set(original_index) & pericyte_aliases)
    if "Pericyte" in original_index and present_pericyte_aliases:
        raise ValueError(
            f"scDRS group table contains both final Pericyte and legacy "
            f"Pericyte aliases for {trait}: "
            + ", ".join(present_pericyte_aliases)
        )
    legacy_labels = sorted(set(original_index) & set(legacy_to_final))
    if legacy_labels:
        standardized.index = original_index.map(
            lambda celltype: legacy_to_final.get(celltype, celltype)
        )
        duplicate_labels = sorted(
            set(standardized.index[standardized.index.duplicated(keep=False)])
        )
        if duplicate_labels:
            raise ValueError(
                f"Legacy-to-final Celltype conversion creates duplicate rows for {trait}: "
                + ", ".join(duplicate_labels)
            )
    else:
        standardized.index = original_index

    available_group_celltypes = set(standardized.index)
    missing_group_celltypes = [
        celltype
        for celltype in atlas_celltype_order
        if celltype not in available_group_celltypes
    ]
    unexpected_group_celltypes = sorted(
        available_group_celltypes - set(atlas_celltype_order)
    )
    if missing_group_celltypes or unexpected_group_celltypes:
        details = []
        if missing_group_celltypes:
            details.append("missing: " + ", ".join(missing_group_celltypes))
        if unexpected_group_celltypes:
            details.append("unexpected: " + ", ".join(unexpected_group_celltypes))
        raise ValueError(
            f"scDRS Celltype rows do not match the final atlas interface for {trait} ("
            + "; ".join(details)
            + ")"
        )

    return standardized.loc[atlas_celltype_order]


dict_df_stats = {
    trait: standardize_group_stats_index(df_stats, trait)
    for trait, df_stats in dict_df_stats.items()
}

# ------------ 3. 绘制分组统计热图（适配你的数据） ------------
fig, ax = scdrs.util.plot_group_stats(
    dict_df_stats=dict_df_stats,
    plot_kws={
        "vmax": 0.2,
        "cb_fraction":0.12
    }
)

fig.savefig("basic_group_stats.png", dpi=300,bbox_inches='tight')
plt.close(fig)

##########################
target_celltypes = [
    "OFB-SMC",
    "Fibroblast",
    "non-OFB-SMC",
]
available_celltypes = set(adata.obs["Celltype"].dropna().astype(str))
missing_target_celltypes = [
    celltype
    for celltype in target_celltypes
    if celltype not in available_celltypes
]

if missing_target_celltypes:
    raise ValueError(
        "Required scDRS cell types are missing: "
        + ", ".join(missing_target_celltypes)
    )

adata_ca1 = adata[adata.obs["Celltype"].isin(target_celltypes)].copy()
adata_ca1.obs["Celltype"] = pd.Categorical(
    adata_ca1.obs["Celltype"],
    categories=target_celltypes,
    ordered=True,
)

for trait in dict_score:
    adata_ca1.obs[trait] = dict_score[trait]["norm_score"]

# 绘制UMAP图，颜色映射为疾病评分
sc.pl.umap(
    adata_ca1,
    color=dict_score.keys(),  # 按疾病评分着色
    color_map="RdBu_r",       # 红蓝渐变（负值蓝，正值红）
    vmin=-5, vmax=5,          # 颜色范围限制
    ncols=3, s=20,            # 排版（3列，点大小20）
    show=False,
)
plt.savefig(
    str(
        configured_path(
            PATHS,
            "SMC_RESULTS_DIR",
            "scdrs",
            "focused_celltypes_umap_plot.png",
        )
    ),
    dpi=300,
    bbox_inches="tight",
)
plt.close()
