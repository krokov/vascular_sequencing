# run_pipeline.py
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from src.plotting import load_and_align_data, compute_color_bounds, draw_group_projection

def main():
    parser = argparse.ArgumentParser(description="ESAI UMAP Overlay Pipeline")
    parser.add_argument("--umap-csv", type=str, required=True, help="Path to cell-type UMAP coordinates")
    parser.add_argument("--esai-csv", type=str, required=True, help="Path to computed ESAI value tables")
    parser.add_argument("--outdir", type=str, default="outputs/plots", help="Target output directory for rendered figures")
    parser.add_argument("--multi-page", action="store_true", help="Compile a consolidated multi-page PDF document")
    args = parser.parse_args()

    umap_path = Path(args.umap_csv)
    esai_path = Path(args.esai_csv)
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Processing and aligning dataset metrics...")
    merged_df, group_col = load_and_align_data(umap_path, esai_path)
    
    unique_groups = sorted(merged_df["group"].dropna().unique())
    vmin, vmax = compute_color_bounds(merged_df["ESAI"].values)
    print(f"[INFO] Determined shared visualization scale: vmin={vmin:.4f}, vmax={vmax:.4f}")

    if args.multi_page:
        pdf_target = output_dir / "ESAI_UMAP_Combined_Report.pdf"
        print(f"[INFO] Rendering consolidated multi-page document to: {pdf_target}")
        with PdfPages(pdf_target) as pdf:
            for group in unique_groups:
                group_slice = merged_df[merged_df["group"] == group]
                if group_slice.empty:
                    continue
                    
                fig, ax = plt.subplots(figsize=(6.5, 6))
                draw_group_projection(ax, merged_df, group_slice, vmin, vmax, title=f"ESAI Overlay — {group}")
                pdf.savefig(fig, dpi=300, bbox_inches="tight")
                plt.close(fig)

    print(f"[INFO] Dispatching individual high-res PNG plots to: {output_dir}")

    for group in unique_groups:
        group_slice = merged_df[merged_df["group"] == group]
        if group_slice.empty:
            continue
            
        fig, ax = plt.subplots(figsize=(6.5, 6))
        draw_group_projection(ax, merged_df, group_slice, vmin, vmax, title=f"ESAI Profile — {group}")
        
        #file names
        clean_name = str(group).strip().lower().replace(" ", "_")
        fig.savefig(output_dir / f"umap_{clean_name}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("[SUCCESS] Pipeline execution cycle successfully finalized.")

if __name__ == "__main__":
    main()