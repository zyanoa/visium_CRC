#!/usr/bin/env python3
# ==============================================================================
# Script: 03_run_cell2location.py
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Load Visium spatial transcriptomics data
#   2. Load scRNA-seq reference data
#   3. Train cell2location reference regression model
#   4. Map reference cell types onto Visium spots
#   5. Export cell abundance matrix and spatial plots
# ==============================================================================

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

rcParams["pdf.fonttype"] = 42


def parse_args():
    parser = argparse.ArgumentParser(description="Run cell2location on Visium data")
    parser.add_argument("--input_sc", type=str, required=True, help="Path to scRNA-seq reference data")
    parser.add_argument("--sc_type", type=str, choices=["csv", "10x"], required=True, help="Reference input type")
    parser.add_argument("--input_sp", type=str, required=True, help="Path to Visium spatial data")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--sample", type=str, required=True, help="Sample name")
    parser.add_argument("--reference_samples", type=str, default="C112,C124,C139,C155,C157,C162,C163,C165,C170,C173", help="Comma-separated CRC reference sample IDs to retain")
    parser.add_argument("--tumor_type_label", type=str, default="T", help="Value in adata_ref.obs['type'] indicating tumor samples; T means Tumor in the GSE178341 reference annotation")
    parser.add_argument("--detection_alpha", type=float, default=20, help="cell2location detection_alpha parameter")
    return parser.parse_args()


def load_spatial_data(spatial_dir: str, sample_name: str):
    adata_vis = sc.read_visium(spatial_dir)
    adata_vis.obs["sample"] = sample_name

    adata_vis.var["SYMBOL"] = adata_vis.var_names
    adata_vis.var.set_index("gene_ids", drop=True, inplace=True)

    adata_vis.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values].copy()

    return adata_vis


def build_symbol_to_ensembl_map(adata_vis):
    return {
        adata_vis.var["SYMBOL"][i]: adata_vis.var.index[i]
        for i in range(adata_vis.var.shape[0])
    }


def load_reference_data(scrna_path: str, sc_type: str, symbol_to_ensg: dict):
    if sc_type == "10x":
        adata_ref = sc.read(scrna_path)
        adata_ref.var["SYMBOL"] = adata_ref.var_names
        adata_ref.var["gene_ids"] = (
            adata_ref.var["gene_ids"]
            .str.split(".").str[0]
            .str.split("_").str[0]
        )
        adata_ref.var.set_index("gene_ids", drop=True, inplace=True)
        adata_ref = adata_ref[:, adata_ref.var.index.str.match("^ENSG")].copy()
        adata_ref.var_names_make_unique()
    else:
        adata = pd.read_csv(scrna_path, index_col=0).T
        tmp_csv = "tmp_reference.csv"
        adata.to_csv(tmp_csv)
        adata_ref = sc.read_csv(tmp_csv)
        os.remove(tmp_csv)

        shared_symbols = [g for g in adata_ref.var.index if g in symbol_to_ensg.keys()]
        adata_ref = adata_ref[:, shared_symbols].copy()
        adata_ref.var["gene_ids"] = [symbol_to_ensg[g] for g in adata_ref.var.index]
        adata_ref.var["SYMBOL"] = adata_ref.var.index
        adata_ref.var.index = adata_ref.var["gene_ids"]

    return adata_ref


def subset_reference_cells(adata_ref, reference_samples, tumor_type_label):
    """Subset the scRNA-seq reference used for the manuscript workflow.

    Important note:
    In this reference annotation, adata_ref.obs["type"] == "T" denotes Tumor,
    not T cells. Cell-type labels for deconvolution are stored in clMidwayPr.
    """
    if "orig.ident" not in adata_ref.obs.columns:
        raise KeyError("adata_ref.obs must contain 'orig.ident'.")
    if "type" not in adata_ref.obs.columns:
        raise KeyError("adata_ref.obs must contain 'type', where T denotes Tumor.")

    sample_ids = [x.strip() for x in reference_samples.split(",") if x.strip()]
    mask = (
        adata_ref.obs["orig.ident"].isin(sample_ids)
    ) & (adata_ref.obs["type"] == tumor_type_label)

    print(f"[INFO] Retained {int(mask.sum())} reference cells from tumor samples: {sample_ids}")
    return adata_ref[mask].copy()


def train_reference_model(adata_ref):
    selected = filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12
    )
    adata_ref = adata_ref[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key="orig.ident",
        labels_key="clMidwayPr"
    )

    model = RegressionModel(adata_ref)
    model.train(max_epochs=1000)

    adata_ref = model.export_posterior(adata_ref)
    adata_ref = model.export_posterior(
        adata_ref,
        use_quantiles=True,
        add_to_varm=["q05", "q50", "q95", "q0001"]
    )

    return adata_ref, model


def extract_reference_signatures(adata_ref):
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = adata_ref.var[
            [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
        ].copy()

    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
    return inf_aver


def run_cell2location_mapping(adata_vis, inf_aver, detection_alpha):
    shared_genes = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

    model = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=20,
        detection_alpha=detection_alpha
    )
    model.train(max_epochs=30000, batch_size=None, train_size=1)

    adata_vis = model.export_posterior(adata_vis)
    return adata_vis, model


def save_outputs(adata_vis, adata_ref, outdir, sample_name):
    factor_names = list(adata_vis.uns["mod"]["factor_names"])
    adata_vis.obs[factor_names] = adata_vis.obsm["q05_cell_abundance_w_sf"]

    abundance_file = os.path.join(outdir, f"{sample_name}.spatial.deconvolution.csv")
    adata_vis.obs[factor_names].to_csv(abundance_file)

    ref_file = os.path.join(outdir, f"{sample_name}.reference_model.h5ad")
    adata_ref.write(ref_file)

    sc.pl.spatial(
        adata_vis,
        cmap="bwr",
        color=factor_names,
        ncols=4,
        size=1.3,
        img_key="hires",
        vmin=0,
        vmax="p99.2",
        show=False
    )
    plt.savefig(
        os.path.join(outdir, f"{sample_name}.spatial.celltype.png"),
        bbox_inches="tight",
        dpi=300
    )
    plt.close()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=300, facecolor="white")

    adata_vis = load_spatial_data(args.input_sp, args.sample)
    symbol_to_ensg = build_symbol_to_ensembl_map(adata_vis)

    adata_ref = load_reference_data(args.input_sc, args.sc_type, symbol_to_ensg)
    adata_ref = subset_reference_cells(adata_ref, args.reference_samples, args.tumor_type_label)

    adata_ref, ref_model = train_reference_model(adata_ref)
    inf_aver = extract_reference_signatures(adata_ref)

    adata_vis, spatial_model = run_cell2location_mapping(adata_vis, inf_aver, args.detection_alpha)
    save_outputs(adata_vis, adata_ref, args.outdir, args.sample)


if __name__ == "__main__":
    main()