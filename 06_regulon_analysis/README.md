# SCENIC regulon analysis for Stroma1_IGLC1 and Stroma4_IGHM

This directory contains the manuscript-facing pySCENIC workflow used to examine
transcription-factor regulon activity in the two humoral/TLS-associated spatial
compartments, **Stroma1_IGLC1** and **Stroma4_IGHM**, across dMMR and pMMR CRC.

## Analysis design

The pySCENIC input contains Visium spots from four groups:

- `dMMR_Stroma1_IGLC1`
- `pMMR_Stroma1_IGLC1`
- `dMMR_Stroma4_IGHM`
- `pMMR_Stroma4_IGHM`

The workflow uses raw counts from the `Spatial` assay. To match the analysis
used during manuscript development, genes are retained when they are detected
in at least `max(20, 1% of selected spots)` and have at least 50 total counts
across the selected spots.

The core pySCENIC workflow is:

1. GRNBoost2 network inference.
2. cisTarget motif enrichment/pruning.
3. AUCell regulon activity scoring.
4. Regulon specificity score (RSS) calculation across the four MMR/region groups.

Manuscript-focused outputs highlight:

- **Stroma4_IGHM / B-cell and germinal-center regulons:** `BCL6`, `IRF8`, `PAX5`, `SPIB`.
- **Stroma1_IGLC1 / plasma-cell regulons:** `IRF4`, `XBP1`, `PRDM1`.

The code reports missing target regulons rather than substituting unrelated
transcription factors.

## Scripts

- `01_prepare_Stroma1_4_expr.R` — extract Stroma1/Stroma4 raw-count matrix and metadata.
- `02_csv_to_loom.py` — convert the historical gene-by-spot CSV input to loom.
- `03_run_pySCENIC_Stroma1_4.sh` — run GRNBoost2, cisTarget, and AUCell.
- `04_analyze_SCENIC_RSS_Stroma1_4.R` — calculate four-group RSS and generate Fig. 3f/3h-style ranking plots.
- `05_plot_SCENIC_spatial_Stroma1_4.R` — map target regulon AUCell scores back to Visium sections without overwriting the Spatial assay.

## Required pySCENIC resources

The workflow expects the following resources to be provided locally (paths are configurable with environment variables):

- `hs_hgnc_tfs.txt`
- `hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`

The manuscript analysis used pySCENIC/SCENIC with AUCell. The exact package and
environment versions should be recorded from the generated session/version logs.

## Example

From the repository root:

```bash
export CRC_MMR_ST_PROJECT="$PWD"
export SCENIC_SEURAT_RDS="$PWD/data/processed/st_obj.rds"
export SCENIC_TF_LIST="/path/to/hs_hgnc_tfs.txt"
export SCENIC_RANKING_DB="/path/to/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
export SCENIC_MOTIF_ANNOT="/path/to/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

sbatch 06_regulon_analysis/03_run_pySCENIC_Stroma1_4.sh
```

After pySCENIC finishes:

```bash
WORK="$PWD/results/06_regulon_analysis/Stroma1_4"

Rscript --vanilla 06_regulon_analysis/04_analyze_SCENIC_RSS_Stroma1_4.R \
  --loom "$WORK/output/sce_SCENIC.loom" \
  --metadata "$WORK/input/Stroma1_4_metadata.csv" \
  --outdir "$WORK/output"

Rscript --vanilla 06_regulon_analysis/05_plot_SCENIC_spatial_Stroma1_4.R \
  --seurat_rds "$PWD/data/processed/st_obj.rds" \
  --auc_csv "$WORK/output/SCENIC_regulon_AUC_Stroma1_4.csv" \
  --outdir "$WORK/output/spatial"
```

The spatial plotting script defaults to the representative sections used in the
historical SCENIC figure workflow (`PT9_2`, `PT9_3`, `PT34_1`, `PT55_2`) but the
set can be changed with `--sections`.
