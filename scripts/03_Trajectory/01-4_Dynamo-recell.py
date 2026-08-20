# =============================================================================
# Title: 01-4_Dynamo-recell.py
# Purpose: Figure 4 and supplementary
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `Merge.h5ad`; processed SMC Dynamo h5ad; metadata/UMAP tables
# Generated outputs: Velocity, streamline, topology and vector-field plots
# Upstream dependencies: `Merge.h5ad`; processed SMC Dynamo h5ad; metadata/UMAP tables
# Downstream consumers: Figure 4 and supplementary
# Configuration keys: DATA_PROCESSED_DIR, TRAJECTORY_RESULTS_DIR, MONOCLE2_PRIVATE_RDS, MONOCLE2_PUBLIC_RDS
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/03_Trajectory/01-4_Dynamo-recell.py
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_shared = _ReleasePath(__file__).resolve().parents[1] / "_shared"
_release_sys.path.insert(0, str(_release_shared))
from paths import configured_path, load_project_paths
PATHS = load_project_paths()
# -*- coding: utf-8 -*-
import os
import sys
import re
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import random
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")   # 鏈嶅姟鍣ㄦ棤鍥惧舰鐣岄潰鏃朵娇鐢?Agg锛岄伩鍏?X11/TkAgg 鎶ラ敊
import matplotlib.pyplot as plt

import scanpy as sc
import dynamo as dyn
import seaborn as sns
from dynamo.preprocessing import Preprocessor

def set_global_seed(seed=42):
    """设置全局随机种子，确保结果可复现"""
    import random
    import numpy as np
    import os
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    # 可选：PyTorch
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass
    # 可选：TensorFlow
    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
    except ImportError:
        pass


set_global_seed(seed=42)
# =========================
# Basic settings
# =========================

model = "monocle"

outpath = str(configured_path(PATHS, "TRAJECTORY_RESULTS_DIR", "dynamo", "monocle"))
readpath = str(configured_path(PATHS, "TRAJECTORY_RESULTS_DIR", "dynamo"))

os.makedirs(outpath, exist_ok=True)
os.chdir(readpath)

celltype_key = "Celltype"          
plot_celltype_key = "Celltype_plot"  


# =========================
# Original color map for analysis labels
# =========================

pancreas_cluster_cmap = {
    "Contractile SMC": "#3361A5",
    "Transitional SMC": "#248AF3",
    "Macrophage_like_SMC_1": "#FDBF6F",
    "Macrophage_like_SMC_2": "#FF7F00",
    "Fibromyocyte-like SMC": "#A6CEE3",
    "OFB-SMC": "#B81136",
}


# =========================
# New plot labels and color map
# =========================

rename_map = {
    "Contractile SMC": "Contractile SMC",
    "Transitional SMC": "Transitional SMC",
    "Macrophage_like_SMC_1": "Macrophage-like SMC-1",
    "Macrophage_like_SMC_2": "Macrophage-like SMC-2",
    "Fibromyocyte-like SMC": "Fibromyocyte-like SMC",
    "OFB-SMC": "OFB-SMC",
}

plot_cluster_cmap = {
    "Contractile SMC": "#add488",
    "Transitional SMC": "#add488",
    "Macrophage-like SMC-1": "#f5ba70",
    "Macrophage-like SMC-2": "#ef7c1a",
    "Fibromyocyte-like SMC": "#a3c9dc",
    "OFB-SMC": "#389937",
}


# =========================
# Read data
# =========================

adata = sc.read_h5ad(
    str(configured_path(PATHS, "DYNAMO_DATA_DIR", "Merge.h5ad"))
)
sub = sc.read_h5ad(str(configured_path(PATHS, "TRAJECTORY_RESULTS_DIR", "dynamo", "monocle", "Dy_smc-monocle_pp.h5ad")))
#鍓嶆湡宸茬粡鍋氬畬浜嗗垎鏋愪簡锛屾墍浠ョ洿鎺ョ敤
meta = pd.read_csv(
    os.path.join(readpath, "sub-meta-dynamo.xls"),
    index_col=0,
    sep="\t"
)

adata.obs.index = adata.obs.index.str.replace(":", "_", regex=False)

sub = adata[meta.index].copy()
sub.obs = meta.copy()

umap = pd.read_csv(
    os.path.join(readpath, "sub-umap-dynamo.xls"),
    sep="\t",
    index_col=0
)

sub.obs[["UMAP_1", "UMAP_2"]] = umap[["UMAP_1", "UMAP_2"]]
sub.obsm["X_umap"] = sub.obs[["UMAP_1", "UMAP_2"]].values


# =========================
# Check and preprocess
# =========================

if not sub.var_names.is_unique:
    duplicated_var_names = sub.var_names[sub.var_names.duplicated()].unique()
    print(f"{duplicated_var_names}")

sub.var_names_make_unique(join="-")

sc.pp.calculate_qc_metrics(
    sub,
    expr_type="unspliced",
    layer="unspliced",
    inplace=True
)

sc.pp.calculate_qc_metrics(
    sub,
    expr_type="spliced",
    layer="spliced",
    inplace=True
)

sc.pl.violin(sub, ["n_genes_by_spliced", "n_genes_by_unspliced"])
plt.savefig(os.path.join(outpath, "unspliced_spliced_violin.n_genes.png"), dpi=300, bbox_inches="tight")
plt.close()

sc.pl.violin(sub, ["log1p_total_spliced", "log1p_total_unspliced"])
plt.savefig(os.path.join(outpath, "unspliced_spliced_violin.log1p_total.png"), dpi=300, bbox_inches="tight")
plt.close()

sc.pl.violin(sub, ["total_spliced", "total_unspliced"])
plt.savefig(os.path.join(outpath, "unspliced_splice_violin.total.png"), dpi=300, bbox_inches="tight")
plt.close()


# =========================
# Dynamo preprocessing
# =========================

preprocessor = Preprocessor()
preprocessor.preprocess_adata(sub, recipe=model)


# =========================
# Keep original Celltype for analysis
# Create Celltype_plot only for plotting
# =========================

sub.obs[celltype_key] = sub.obs[celltype_key].astype(str)
sub.obs[celltype_key] = sub.obs[celltype_key].str.replace("-", "_", regex=False)
sub.obs[celltype_key] = sub.obs[celltype_key].str.replace(" ", "_", regex=False)

sub.obs[plot_celltype_key] = sub.obs[celltype_key].replace(rename_map)

print("Original Celltype categories:")
print(sub.obs[celltype_key].value_counts())

print("\nPlot Celltype categories:")
print(sub.obs[plot_celltype_key].value_counts())

missing_plot_colors = set(sub.obs[plot_celltype_key].unique()) - set(plot_cluster_cmap.keys())
if len(missing_plot_colors) > 0:
    print("\nWarning: these plot labels do not have colors in plot_cluster_cmap:")
    print(missing_plot_colors)


# =========================
# Dynamics
# =========================

dyn.tl.dynamics(sub, model="stochastic")

dyn.tl.reduceDimension(sub, basis="pca")


# =========================
# UMAP with new labels
# =========================

dyn.pl.umap(
    sub,
    color=plot_celltype_key,
    pointsize=0.2,
    alpha=0.3,
    color_key=plot_cluster_cmap,
    show_legend="False",
    inset_dict={
        "width": "0.05",
        "height": "0.25",
        "loc": "lower right",
        "bbox_to_anchor": (0.8, 0.1),
    },
    show_arrowed_spines=True,
    theme="darkred",
    save_show_or_return="show",
)


# =========================
# Confidence and velocities
# Important: keep Celltype here, not Celltype_plot
# =========================

dyn.tl.gene_wise_confidence(
    sub,
    group=celltype_key,
    lineage_dict={"Contractile SMC": ["OFB-SMC"]},
)

dyn.tl.cell_velocities(
    sub,
    method="fp",
    basis="umap",
    enforce=True,
    transition_genes=list(sub.var_names[sub.var.use_for_pca]),
)

dyn.tl.cell_wise_confidence(sub)

dyn.tl.confident_cell_velocities(
    sub,
    group=celltype_key,
    lineage_dict={"Contractile SMC": ["OFB-SMC"]},
    only_transition_genes=True,
    basis="umap",
)

dyn.tl.cell_wise_confidence(sub)

dyn.tl.cell_velocities(sub, basis="pca")

dyn.vf.VectorField(sub, basis="pca", calc_jacobian=True)
dyn.vf.VectorField(sub, basis="umap", M=1000, pot_curl_div=True)


# =========================
# Streamline plot with new labels
# =========================

dyn.pl.streamline_plot(
    sub,
    color=[plot_celltype_key],
    basis="umap",
    color_key=plot_cluster_cmap,
    show_legend="False",
    save_show_or_return="show",
)
dyn.pl.streamline_plot(sub, color=[plot_celltype_key], basis="umap",color_key = plot_cluster_cmap,show_legend='False',save_show_or_return="save", adjust_text=True,text_force_strength=0.5,text_min_distance=5,text_fontsize=6,save_kwargs={"path":outpath,"prefix": f"streamline_plotff","dpi": 900,"ext": 'svg', "transparent": True, "close": True, "verbose": True})


# =========================
# Velocity / acceleration / curvature ranking
# Important: keep Celltype here, not Celltype_plot
# =========================

dyn.vf.rank_velocity_genes(
    sub,
    groups=celltype_key,
    vkey="velocity_S",
)

rank_speed = sub.uns["rank_velocity_S"]
rank_abs_speed = sub.uns["rank_abs_velocity_S"]

dyn.vf.acceleration(sub, basis="pca")

dyn.vf.rank_acceleration_genes(
    sub,
    groups=celltype_key,
    akey="acceleration",
    prefix_store="rank",
)

rank_acceleration = sub.uns["rank_acceleration"]
rank_abs_acceleration = sub.uns["rank_abs_acceleration"]

dyn.vf.curvature(sub, basis="pca")

dyn.vf.rank_curvature_genes(
    sub,
    groups=celltype_key,
)

sub.obs["speed_pca"] = np.linalg.norm(sub.obsm["velocity_pca"], axis=1)

dyn.vf.divergence(sub, basis="pca")


# =========================
# Topography with new labels
# =========================

dyn.pl.topography(
    sub,
    basis="umap",
    background="white",
    color=[plot_celltype_key],
    streamline_color="black",
    color_key=plot_cluster_cmap,
    save_show_or_return="show",
)


# =========================
# Potential plot
# =========================

dyn.pl.umap(
    sub,
    color="umap_ddhodge_potential",
    frontier=True,
    save_show_or_return="show",
)


# =========================
# Gene expression along potential with new labels
# =========================

fig3_si5 = ["CREB3L2", "RUNX2", "BMP2", "SPP1", "GATA6", "SIRT6", "NR4A3"]

genes = [
    "MYO1D",
    "PKD1",
    "SMAD7",
    "PLXDC2",
    "ANXA4",
    "AIG1",
    "PGRMC2",
    "SLC17A5",
    "TMEM50B",
]

for gene in genes:
    dyn.pl.scatters(
        sub,
        x="umap_ddhodge_potential",
        y=gene,
        layer="X_spliced",
        color=plot_celltype_key,
        pointsize=0.25,
        alpha=0.8,
        background="white",
        color_key=plot_cluster_cmap,
        show_legend="False",
        save_show_or_return="save",
        save_kwargs={
            "path": outpath,
            "prefix": f"{gene}_pl_scatters",
            "dpi": None,
            "ext": "png",
            "transparent": True,
            "close": True,
            "verbose": True,
        },
    )


# =========================
# Final velocity gene ranking
# Important: keep Celltype here, not Celltype_plot
# =========================

dyn.vf.rank_velocity_genes(
    sub,
    groups=celltype_key,
    vkey="velocity_S",
)

# Save processed object if needed
# sub.write(os.path.join(outpath, "Dy_sub_monocle_pp_renamed_plot_labels.h5ad"))
