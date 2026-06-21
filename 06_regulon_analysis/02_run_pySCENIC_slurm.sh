#!/bin/bash
#SBATCH --job-name=pySCENIC_Stroma3_4
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --error=%j.err
#SBATCH --output=%j.out

set -euo pipefail

# Usage:
#   sbatch 06_regulon_analysis/02_run_pySCENIC_slurm.sh \
#     /path/to/scenic_workdir \
#     /path/to/hs_hgnc_tfs.txt \
#     /path/to/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
#     /path/to/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

WORKDIR=${1:-"results/06_regulon_analysis/SCENIC_Stroma3_4"}
TF_LIST=${2:-"data/external/scenic/hs_hgnc_tfs.txt"}
MOTIF_ANNOT=${3:-"data/external/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"}
RANKING_DB=${4:-"data/external/scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"}

source ~/.bash_profile
conda activate pyscenic

cd "${WORKDIR}"

pyscenic grn --num_workers 10 \
  --sparse \
  --method grnboost2 \
  --output sce.adj.csv \
  sce.loom \
  "${TF_LIST}"

pyscenic ctx --num_workers 10 \
  --output sce.regulons.csv \
  --expression_mtx_fname sce.loom \
  --all_modules \
  --mask_dropouts \
  --mode "dask_multiprocessing" \
  --min_genes 20 \
  --annotations_fname "${MOTIF_ANNOT}" \
  sce.adj.csv \
  "${RANKING_DB}"

pyscenic aucell --num_workers 10 \
  --output sce_SCENIC.loom \
  sce.loom \
  sce.regulons.csv
