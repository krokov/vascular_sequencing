# src/pipeline.py
from pathlib import Path
from typing import Set, List
import pandas as pd

def parse_curated_gene_list(file_path: Path) -> Set[str]:
    """
    Reads a curated list of target marker genes (TXT, CSV, or TSV).
    Maintains formatting and returns a unique set for quick intersections.
    """
    if not file_path.exists():
        print(f"[WARNING] Requested gene file target not found: {file_path}")
        return set()
        
    ext = file_path.suffix.lower()
    try:
        if ext == '.csv':
            df = pd.read_csv(file_path, header=None)
            genes = df[0].dropna().astype(str).tolist()
        elif ext in ['.tsv', '.txt']:
            df = pd.read_csv(file_path, sep='\t', header=None)
            genes = df[0].dropna().astype(str).tolist()
        else:
            #fallback plaintext
            with open(file_path, 'r', encoding='utf-8') as f:
                genes = f.read().splitlines()
                
        #clean
        cleaned_genes = {g.strip() for g in genes if g.strip()}
        print(f"[INFO] Successfully loaded {len(cleaned_genes)} curated target genes.")
        return cleaned_genes
        
    except Exception as e:
        print(f"[ERROR] Failed to extract custom gene file metrics: {e}")
        return set()

def sync_pipeline_metadata(output_dir: Path, sample_names: List[str]) -> pd.DataFrame:
    """
    Generates a tracking matrix framework to verify output consistency
    across all processed experimental runs.
    """
    manifest_path = output_dir / "pipeline_run_manifest.csv"
    records = []
    
    for sample in sample_names:
        records.append({
            "sample_id": sample,
            "data_aligned": "Pending",
            "plots_generated": "False"
        })
        
    df = pd.DataFrame(records)
    df.to_csv(manifest_path, index=False)
    return df