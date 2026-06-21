#!/usr/bin/env python3
"""
Prepare a loom file for pySCENIC from a CSV expression matrix.

Expected input:
  - A gene-by-spot/cell or spot/cell-by-gene CSV expression matrix.
  - By default, rows are observations and columns are genes, matching Scanpy's read_csv behavior.

This script is used for the Stroma3_APOE / Stroma4_IGHM regulon analysis.
"""

import argparse
import os
import numpy as np
import scanpy as sc
import loompy as lp


def parse_args():
    parser = argparse.ArgumentParser(description="Create pySCENIC loom input from CSV expression matrix.")
    parser.add_argument("--input_csv", required=True, help="Input expression CSV, e.g. Stroma3_4.csv")
    parser.add_argument("--output_loom", required=True, help="Output loom file, e.g. sce.loom")
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_csv(args.input_csv)

    row_attrs = {"Gene": np.array(adata.var_names)}
    col_attrs = {"CellID": np.array(adata.obs_names)}

    os.makedirs(os.path.dirname(os.path.abspath(args.output_loom)), exist_ok=True)
    lp.create(args.output_loom, adata.X.transpose(), row_attrs, col_attrs)
    print(f"Wrote loom file: {args.output_loom}")


if __name__ == "__main__":
    main()
