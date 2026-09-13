#!/usr/bin/env python3

"""Convert a gene-by-spot CSV expression matrix to loom format for pySCENIC.

The first column must contain unique gene names and the remaining columns must
contain unique Visium spot IDs. The matrix is kept in gene x spot orientation,
which is the orientation expected by the pySCENIC CLI.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import loompy
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert gene-by-spot CSV expression matrix to loom format for pySCENIC."
    )
    parser.add_argument(
        "--expr",
        required=True,
        type=Path,
        help="Input CSV. First column = gene, remaining columns = spots.",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output loom file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.expr.exists():
        raise FileNotFoundError(f"Expression CSV not found: {args.expr}")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading expression CSV: {args.expr}", flush=True)
    expr = pd.read_csv(args.expr)

    if expr.shape[1] < 2:
        raise ValueError("Expression CSV must contain one gene column and at least one spot column.")

    gene_col = expr.columns[0]
    genes = expr.iloc[:, 0].astype(str).to_numpy()
    mat_df = expr.iloc[:, 1:]
    spots = mat_df.columns.astype(str).to_numpy()

    if pd.Index(genes).duplicated().any():
        dup = pd.Index(genes)[pd.Index(genes).duplicated()].unique().tolist()[:10]
        raise ValueError(f"Duplicated gene names found. Examples: {dup}")

    if pd.Index(spots).duplicated().any():
        dup = pd.Index(spots)[pd.Index(spots).duplicated()].unique().tolist()[:10]
        raise ValueError(f"Duplicated spot IDs found. Examples: {dup}")

    matrix = mat_df.to_numpy(dtype=np.float32, copy=True)

    if not np.isfinite(matrix).all():
        raise ValueError("Expression matrix contains NaN or Inf values.")
    if (matrix < 0).any():
        raise ValueError("Expression matrix contains negative values; raw non-negative counts are expected.")

    print(f"[INFO] Gene column: {gene_col}", flush=True)
    print(f"[INFO] Matrix: {matrix.shape[0]} genes x {matrix.shape[1]} spots", flush=True)

    row_attrs = {"Gene": genes}
    col_attrs = {"CellID": spots, "cell_id": spots}

    if args.out.exists():
        args.out.unlink()

    loompy.create(
        filename=str(args.out),
        layers=matrix,
        row_attrs=row_attrs,
        col_attrs=col_attrs,
    )

    with loompy.connect(str(args.out), mode="r") as ds:
        if ds.shape != matrix.shape:
            raise RuntimeError(f"Loom shape mismatch: expected {matrix.shape}, observed {ds.shape}")
        print(f"[INFO] Loom created: {args.out}", flush=True)
        print(f"[INFO] Loom shape: {ds.shape}", flush=True)
        print(f"[INFO] Row attrs: {list(ds.ra.keys())}", flush=True)
        print(f"[INFO] Col attrs: {list(ds.ca.keys())}", flush=True)


if __name__ == "__main__":
    main()
