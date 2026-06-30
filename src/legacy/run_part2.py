#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ESAI UMAP (sample overlays) — robust plotter
--------------------------------------------
Generates ESAIumap_sample plots per group/condition on a fixed cell-type UMAP.

Features:
- Validates inputs & columns, with helpful errors
- Shared color scale across groups (auto 99th pct or fixed)
- (A) Multi-page PDF (one page per group)
- (B) Separate PNG/PDF files per group
- Skips groups with no data (warns)
- Works for Balanced or Unbalanced inputs

Example:
python run_part2.py \
  --umap-csv data/celltype_umap_coords.csv \
  --esai-csv data/ESAI_celltype_python_BALANCED.csv \
  --group-col condition \
  --groups "Sham" "Sham + EtOH" "Ligated" "Ligated + EtOH" "Ligated + Binge" \
  --outdir outputs/esaiumap_sample/BALANCED \
  --model-name "Balanced" \
  --multi-page \
  --per-group \
  --vmin 0 --vmax 0.5
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# -----------------------------
# Utilities
# -----------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def require_cols(df, cols, name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"[ERROR] {name} is missing required column(s): {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

def slug(s: str) -> str:
    return (
        str(s)
        .strip()
        .lower()
        .replace("+", "plus")
        .replace("&", "and")
        .replace("/", "_")
        .replace(" ", "_")
    )

def safe_mkdir(path: str):
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def compute_shared_bounds(esai_values: np.ndarray, vmin_arg: float = None, vmax_arg: float = None,
                          default_vmax_hint: float = None):
    """Compute shared vmin/vmax. If provided, respect user args; else use robust 99th percentile."""
    # Filter finite
    vals = esai_values[np.isfinite(esai_values)]
    if vals.size == 0:
        return 0.0, 1.0

    # vmin
    vmin = float(vmin_arg) if vmin_arg is not None else 0.0

    # vmax
    if vmax_arg is not None:
        vmax = float(vmax_arg)
    else:
        q99 = float(np.nanquantile(vals, 0.99))
        vmax = q99 if q99 > 0 else (default_vmax_hint if default_vmax_hint is not None else 1.0)

    # Ensure sane
    if vmax <= vmin:
        vmax = vmin + 1e-9

    return vmin, vmax

def draw_one_group(ax, umap_all, umap_group, vmin, vmax, title):
    # Light grey backdrop of all points (shows the map silhouette)
    ax.scatter(
        umap_all["umap_x"], umap_all["umap_y"],
        s=12, c="#e6e6e6", alpha=0.50, linewidths=0
    )

    # Overlay colored ESAI points for this group
    sc = ax.scatter(
        umap_group["umap_x"], umap_group["umap_y"],
        s=28,
        c=umap_group["ESAI"],
        cmap="viridis",
        vmin=vmin, vmax=vmax,
        linewidths=0
    )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_aspect("equal", "box")

    cbar = plt.colorbar(sc, ax=ax, shrink=0.82)
    cbar.set_label("ESAI (per cell type)")

    # Clean up ticks/frames
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Generate ESAIumap_sample plots per group on a fixed cell-type UMAP."
    )
    ap.add_argument("--umap-csv", required=True,
                    help="CSV with cell-type UMAP coordinates. Required columns: source_celltype, umap_x, umap_y")
    ap.add_argument("--esai-csv", required=True,
                    help="CSV with per-group ESAI values per cell type. Must include group and ESAI columns.")
    ap.add_argument("--group-col", default=None,
                    help="Column that defines groups (e.g., 'condition' or 'sample'). "
                         "If omitted, script will auto-detect 'condition' or 'sample'.")
    ap.add_argument("--group-col-rename", default=None,
                    help="Optional: rename the detected/selected group column to this value in outputs.")
    ap.add_argument("--groups", nargs="+", default=None,
                    help="Explicit list of groups to plot in order. If omitted, uses unique sorted values from group column.")
    ap.add_argument("--outdir", default="outputs/esaiumap_sample",
                    help="Output directory for figures.")
    ap.add_argument("--model-name", default=None,
                    help="Label printed in titles (e.g., 'Balanced' or 'Unbalanced').")
    ap.add_argument("--multi-page", action="store_true",
                    help="Write a multi-page ESAIumap_sample.pdf (one page per group).")
    ap.add_argument("--per-group", action="store_true",
                    help="Write separate PNG/PDF files per group.")
    ap.add_argument("--vmin", type=float, default=None,
                    help="Fixed vmin for color scale. Default: auto=0.")
    ap.add_argument("--vmax", type=float, default=None,
                    help="Fixed vmax for color scale. Default: auto=99th percentile.")
    ap.add_argument("--default-vmax-hint", type=float, default=None,
                    help="If auto vmax resolves to 0, use this fallback (e.g., 0.5 for Balanced, 7.0 for Unbalanced).")
    args = ap.parse_args()

    # Basic I/O checks
    if not os.path.exists(args.umap_csv):
        eprint(f"[ERROR] --umap-csv not found: {args.umap_csv}")
        sys.exit(1)
    if not os.path.exists(args.esai_csv):
        eprint(f"[ERROR] --esai-csv not found: {args.esai_csv}")
        sys.exit(1)

    # Load data
    try:
        umap_df = pd.read_csv(args.umap_csv)
    except Exception as ex:
        eprint(f"[ERROR] Failed to read UMAP CSV: {ex}")
        sys.exit(1)

    try:
        esai_df = pd.read_csv(args.esai_csv)
    except Exception as ex:
        eprint(f"[ERROR] Failed to read ESAI CSV: {ex}")
        sys.exit(1)

    # Validate columns in UMAP
    require_cols(umap_df, ["source_celltype", "umap_x", "umap_y"], "UMAP CSV")

    # Determine group column
    gcol = args.group_col
    if gcol is None:
        # auto-detect common names
        for candidate in ["condition", "sample", "group"]:
            if candidate in esai_df.columns:
                gcol = candidate
                break
        if gcol is None:
            eprint("[ERROR] Could not auto-detect group column; provide --group-col (e.g., 'condition' or 'sample').")
            eprint(f"Available columns in ESAI CSV: {list(esai_df.columns)}")
            sys.exit(1)

    # Validate columns in ESAI table
    # Accept ESAI column naming variants
    esai_col_candidates = ["ESAI", "ESAI_c", "esai_value"]
    esai_col = None
    for c in esai_col_candidates:
        if c in esai_df.columns:
            esai_col = c
            break
    if esai_col is None:
        eprint(f"[ERROR] ESAI column not found. Expected one of {esai_col_candidates}. "
               f"Available: {list(esai_df.columns)}")
        sys.exit(1)

    require_cols(esai_df, ["source_celltype", gcol, esai_col], "ESAI CSV")

    # Optional rename of group column
    if args.group_col_rename:
        esai_df = esai_df.rename(columns={gcol: args.group_col_rename})
        gcol = args.group_col_rename

    # Merge: bring UMAP coords onto ESAI rows
    merged = esai_df.merge(umap_df, how="left", on="source_celltype")
    # Warn on unmatched cell types
    n_na = merged["umap_x"].isna().sum()
    if n_na > 0:
        missing_ct = merged.loc[merged["umap_x"].isna(), "source_celltype"].unique().tolist()
        eprint(f"[WARN] {n_na} rows have no UMAP coords (will be dropped). "
               f"Missing cell types (first 10): {missing_ct[:10]}")

    merged = merged.dropna(subset=["umap_x", "umap_y"])
    merged = merged.rename(columns={esai_col: "ESAI"})

    # Determine groups to plot
    if args.groups:
        groups = list(args.groups)
    else:
        groups = sorted(merged[gcol].dropna().unique().tolist())
    if len(groups) == 0:
        eprint("[ERROR] No groups found to plot. Check --group-col and ESAI CSV contents.")
        sys.exit(1)

    # Shared color bounds
    vmin, vmax = compute_shared_bounds(merged["ESAI"].values, args.vmin, args.vmax, args.default_vmax_hint)
    eprint(f"[INFO] Using shared color scale: vmin={vmin:.4g}, vmax={vmax:.4g}")

    # Output dirs
    safe_mkdir(args.outdir)
    pergroup_dir = os.path.join(args.outdir, "per_group")
    if args.per_group:
        safe_mkdir(pergroup_dir)

    # (A) Multi-page PDF
    if args.multi_page:
        mp_path = os.path.join(args.outdir, "ESAIumap_sample.pdf")
        with PdfPages(mp_path) as pdf:
            for grp in groups:
                g = merged.loc[merged[gcol] == grp].copy()
                if g.empty:
                    eprint(f"[WARN] No data for group '{grp}'. Skipping page.")
                    continue

                fig, ax = plt.subplots(figsize=(6.5, 6))
                title = f"ESAIumap_sample — {grp}" + (f" [{args.model_name}]" if args.model_name else "")
                draw_one_group(ax, merged, g, vmin, vmax, title)
                pdf.savefig(fig, dpi=300, bbox_inches="tight")
                plt.close(fig)

        eprint(f"[OK] Wrote multi-page PDF: {mp_path}")

    # (B) Separate files per group
    if args.per_group:
        for grp in groups:
            g = merged.loc[merged[gcol] == grp].copy()
            if g.empty:
                eprint(f"[WARN] No data for group '{grp}'. Skipping files.")
                continue

            fig, ax = plt.subplots(figsize=(6.5, 6))
            title = f"ESAIumap_sample — {grp}" + (f" [{args.model_name}]" if args.model_name else "")
            draw_one_group(ax, merged, g, vmin, vmax, title)

            base = os.path.join(pergroup_dir, f"{slug(grp)}_ESAIumap_sample")
            fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
            fig.savefig(base + ".pdf", dpi=300, bbox_inches="tight")
            plt.close(fig)
            eprint(f"[OK] Wrote: {base}.png / .pdf")

    # If neither option specified, produce at least the multi-page PDF by default
    if not args.multi_page and not args.per_group:
        eprint("[INFO] Neither --multi-page nor --per-group specified. "
               "Defaulting to multi-page PDF output.")
        mp_path = os.path.join(args.outdir, "ESAIumap_sample.pdf")
        with PdfPages(mp_path) as pdf:
            for grp in groups:
                g = merged.loc[merged[gcol] == grp].copy()
                if g.empty:
                    eprint(f"[WARN] No data for group '{grp}'. Skipping page.")
                    continue
                fig, ax = plt.subplots(figsize=(6.5, 6))
                title = f"ESAIumap_sample — {grp}" + (f" [{args.model_name}]" if args.model_name else "")
                draw_one_group(ax, merged, g, vmin, vmax, title)
                pdf.savefig(fig, dpi=300, bbox_inches="tight")
                plt.close(fig)
        eprint(f"[OK] Wrote multi-page PDF: {mp_path}")

    eprint("[DONE] ESAIumap_sample generation complete.")

if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        eprint(f"[FATAL] {ex}")
        sys.exit(1)
