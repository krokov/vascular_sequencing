# src/plotting.py
from pathlib import Path
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages

def load_and_align_data(umap_csv: Path, esai_csv: Path) -> Tuple[pd.DataFrame, str]:
    """Reads UMAP and ESAI data coordinates, auto-detects groupings, and merges them."""
    umap_df = pd.read_csv(umap_csv)
    esai_df = pd.read_csv(esai_csv)
    
    required_umap = {"source_celltype", "umap_x", "umap_y"}
    if not required_umap.issubset(umap_df.columns):
        raise ValueError(f"UMAP file missing required structural coordinates: {required_umap - set(umap_df.columns)}")
    group_candidates = ["condition", "sample", "group"]
    group_col = next((col for col in group_candidates if col in esai_df.columns), None)
    if not group_col:
        raise KeyError(f"Could not isolate grouping column. Expected one of: {group_candidates}")
    value_col = next((col for col in ["ESAI", "ESAI_c", "esai_value"] if col in esai_df.columns), None)
    if not value_col:
        raise KeyError("Could not locate vector score metrics columns in ESAI file.")
        
    merged = esai_df.merge(umap_df, on="source_celltype", how="left").dropna(subset=["umap_x", "umap_y"])
    merged = merged.rename(columns={value_col: "ESAI"})
    
    return merged, group_col

def compute_color_bounds(values: np.ndarray, vmin: float = 0.0, vmax: Optional[float] = None) -> Tuple[float, float]:
    """Establishes deterministic shared upper and lower color ranges."""
    finite_vals = values[np.isfinite(values)]
    if finite_vals.size == 0:
        return vmin, vmin + 1.0
        
    actual_vmax = float(vmax) if vmax is not None else float(np.nanquantile(finite_vals, 0.99))
    #0 divisor prevention
    return vmin, max(actual_vmax, vmin + 1e-9)

def draw_group_projection(ax: Axes, all_data: pd.DataFrame, group_data: pd.DataFrame, vmin: float, vmax: float, title: str):
    """Generates a clean canvas layout overlaying specific experimental cohorts."""
    ax.scatter(all_data["umap_x"], all_data["umap_y"], s=12, c="#e6e6e6", alpha=0.5, linewidths=0)
    
    sc = ax.scatter(group_data["umap_x"], group_data["umap_y"], s=28, c=group_data["ESAI"], 
                    cmap="viridis", vmin=vmin, vmax=vmax, linewidths=0)
    
    ax.set_title(title, fontsize=12)
    ax.set_aspect("equal", "box")
    
    cbar = plt.colorbar(sc, ax=ax, shrink=0.82)
    cbar.set_label("ESAI (per cell type)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])