#!/usr/bin/env python3
"""
COMMOT directional signaling analysis for representative dMMR and pMMR CRC sections.

Purpose
-------
Generate the COMMOT vector-field plots used for the Figure 4E-style comparison of
CCL and CXCL chemokine signaling in representative spatial transcriptomics sections:

    PT18_1 : dMMR
    PT55_2 : pMMR

The script keeps only the analysis required for the manuscript figure:
    - pathways: CCL and CXCL
    - summaries: sender and receiver
    - backgrounds: H&E image and COMMOT summary background

Exploratory pathways and samples from the original notebook were intentionally removed
(COMPLEMENT, WNT, Notch, TGFb, EGF, MK, SPP1, SEMA3, other patient sections, DEG tests).

Expected input
--------------
One AnnData h5ad file per sample under the same COMMOT input directory:

    <input_dir>/PT18_1_adata.h5ad
    <input_dir>/PT55_2_adata.h5ad

The PT18_1 and PT55_2 files should be read from the same regular COMMOT input directory. Do not use SCENIC-specific subdirectories for this analysis.

The files may either already contain COMMOT spatial-communication results, or the script can recompute them using the CellChat secreted signaling database by setting
--run-spatial-communication.

Example
-------
python 05_chemokine_signaling/01_run_COMMOT_CCL_CXCL_representative_samples.py \
    --input-dir data/processed/commot_h5ad \
    --output-dir results/05_chemokine_signaling/COMMOT_Fig4E \
    --run-spatial-communication

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import commot as ct


DEFAULT_SAMPLES: Dict[str, str] = {
    "PT18_1": "dMMR",
    "PT55_2": "pMMR",
}

DEFAULT_PATHWAYS: Tuple[str, ...] = ("CCL", "CXCL")
DEFAULT_SUMMARIES: Tuple[str, ...] = ("sender", "receiver")
DEFAULT_BACKGROUNDS: Tuple[str, ...] = ("image", "summary")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run/plot COMMOT CCL and CXCL directional signaling for PT18_1 and PT55_2."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing <sample>_adata.h5ad files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for COMMOT plots and processed h5ad files.",
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        default=[f"{sample}:{group}" for sample, group in DEFAULT_SAMPLES.items()],
        help="Sample definitions as sample:group. Default: PT18_1:dMMR PT55_2:pMMR.",
    )
    parser.add_argument(
        "--cluster-key",
        default="leiden",
        help="obs column used for background cluster coloring. Default: leiden.",
    )
    parser.add_argument(
        "--run-spatial-communication",
        action="store_true",
        help="Recompute COMMOT spatial communication before plotting. Use this for raw/prepared h5ad inputs.",
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=500.0,
        help="Spatial distance threshold in microns for COMMOT spatial_communication. Default: 500.",
    )
    parser.add_argument(
        "--min-cell-pct",
        type=float,
        default=0.05,
        help="Minimum fraction of spots expressing ligand/receptor for LR filtering. Default: 0.05.",
    )
    parser.add_argument(
        "--k-direction",
        type=int,
        default=5,
        help="k parameter for ct.tl.communication_direction. Default: 5.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=4e-6,
        help="Arrow scaling for ct.pl.plot_cell_communication. Default: 4e-6.",
    )
    parser.add_argument(
        "--grid-density",
        type=float,
        default=0.4,
        help="Grid density for vector-field plotting. Default: 0.4.",
    )
    parser.add_argument(
        "--node-size",
        type=float,
        default=10,
        help="Node size for vector-field plotting. Default: 10.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output image resolution. Default: 600 dpi.",
    )
    parser.add_argument(
        "--save-h5ad",
        action="store_true",
        help="Save h5ad files after COMMOT direction calculation.",
    )
    return parser.parse_args()


def parse_sample_definitions(sample_defs: Iterable[str]) -> Dict[str, str]:
    samples: Dict[str, str] = {}
    for item in sample_defs:
        if ":" not in item:
            raise ValueError(f"Invalid sample definition '{item}'. Use sample:group, e.g. PT18_1:dMMR.")
        sample, group = item.split(":", 1)
        samples[sample] = group
    return samples


def load_ligand_receptor_database() -> pd.DataFrame:
    """Load CellChat secreted signaling ligand-receptor database for COMMOT."""
    return ct.pp.ligand_receptor_database(
        species="human",
        signaling_type="Secreted Signaling",
        database="CellChat",
    )


def maybe_run_spatial_communication(
    adata,
    df_cellchat: pd.DataFrame,
    distance_threshold: float,
    min_cell_pct: float,
):
    """Run COMMOT spatial communication with CCL/CXCL-focused LR database."""
    df_cellchat_filtered = ct.pp.filter_lr_database(
        df_cellchat,
        adata,
        min_cell_pct=min_cell_pct,
    )

    # Column '2' stores pathway names in the COMMOT CellChat LR table.
    df_cellchat_filtered = df_cellchat_filtered[
        df_cellchat_filtered["2"].isin(DEFAULT_PATHWAYS)
    ].copy()

    ct.tl.spatial_communication(
        adata,
        database_name="cellchat",
        df_ligrec=df_cellchat_filtered,
        dis_thr=distance_threshold,
        heteromeric=True,
        pathway_sum=True,
    )
    return df_cellchat_filtered


def plot_one_commot_map(
    adata,
    sample: str,
    group: str,
    pathway: str,
    summary: str,
    background: str,
    output_dir: Path,
    cluster_key: str,
    k_direction: int,
    scale: float,
    grid_density: float,
    node_size: float,
    dpi: int,
):
    """Calculate directionality for one pathway and save one COMMOT vector-field plot."""
    ct.tl.communication_direction(
        adata,
        database_name="cellchat",
        pathway_name=pathway,
        k=k_direction,
    )

    clustering = cluster_key if cluster_key in adata.obs.columns else None
    cmap = "Alphabet" if background == "image" else "Spectral"

    ct.pl.plot_cell_communication(
        adata,
        database_name="cellchat",
        pathway_name=pathway,
        plot_method="grid",
        background_legend=True,
        scale=scale,
        ndsize=node_size,
        grid_density=grid_density,
        summary=summary,
        background=background,
        clustering=clustering,
        cmap=cmap,
        normalize_v=True,
        normalize_v_quantile=0.995,
    )

    output_file = output_dir / f"{sample}_{group}_{pathway}_{summary}_{background}.png"
    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    return output_file


def run_for_sample(
    sample: str,
    group: str,
    input_dir: Path,
    output_dir: Path,
    args: argparse.Namespace,
    df_cellchat: pd.DataFrame,
):
    input_file = input_dir / f"{sample}_adata.h5ad"
    if not input_file.exists():
        raise FileNotFoundError(f"Missing input h5ad: {input_file}")

    print(f"[INFO] Loading {sample} ({group}) from {input_file}")
    adata = sc.read_h5ad(input_file)
    adata.var_names_make_unique()

    if args.run_spatial_communication:
        print(f"[INFO] Running COMMOT spatial_communication for {sample}")
        lr_used = maybe_run_spatial_communication(
            adata=adata,
            df_cellchat=df_cellchat,
            distance_threshold=args.distance_threshold,
            min_cell_pct=args.min_cell_pct,
        )
        lr_file = output_dir / f"{sample}_{group}_CellChat_CCL_CXCL_LR_used.csv"
        lr_used.to_csv(lr_file, index=False)
        print(f"[INFO] Saved filtered LR table: {lr_file}")

    for pathway in DEFAULT_PATHWAYS:
        for summary in DEFAULT_SUMMARIES:
            for background in DEFAULT_BACKGROUNDS:
                print(f"[INFO] Plotting {sample} {pathway} {summary} {background}")
                out = plot_one_commot_map(
                    adata=adata,
                    sample=sample,
                    group=group,
                    pathway=pathway,
                    summary=summary,
                    background=background,
                    output_dir=output_dir,
                    cluster_key=args.cluster_key,
                    k_direction=args.k_direction,
                    scale=args.scale,
                    grid_density=args.grid_density,
                    node_size=args.node_size,
                    dpi=args.dpi,
                )
                print(f"[INFO] Saved {out}")

    if args.save_h5ad:
        h5ad_out = output_dir / f"{sample}_{group}_COMMOT_CCL_CXCL.h5ad"
        adata.write_h5ad(h5ad_out)
        print(f"[INFO] Saved processed h5ad: {h5ad_out}")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    samples = parse_sample_definitions(args.samples)
    df_cellchat = load_ligand_receptor_database()

    for sample, group in samples.items():
        run_for_sample(
            sample=sample,
            group=group,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            args=args,
            df_cellchat=df_cellchat,
        )

    print("[DONE] COMMOT CCL/CXCL representative-sample plotting finished.")


if __name__ == "__main__":
    main()
