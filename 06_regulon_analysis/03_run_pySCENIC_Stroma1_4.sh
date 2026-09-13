#!/usr/bin/env bash
#SBATCH --job-name=scenic_S1_S4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${CRC_MMR_ST_PROJECT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
WORKDIR="${SCENIC_WORKDIR:-${PROJECT_ROOT}/results/06_regulon_analysis/Stroma1_4}"
INPUT_DIR="${WORKDIR}/input"
OUTPUT_DIR="${WORKDIR}/output"
TMP_BASE="${SCENIC_TMPDIR:-${WORKDIR}/tmp}"

SEURAT_RDS="${SCENIC_SEURAT_RDS:-${PROJECT_ROOT}/data/processed/st_obj.rds}"
TF_LIST="${SCENIC_TF_LIST:-${PROJECT_ROOT}/data/external/scenic/hs_hgnc_tfs.txt}"
RANKING_DB="${SCENIC_RANKING_DB:-${PROJECT_ROOT}/data/external/scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"
MOTIF_ANNOT="${SCENIC_MOTIF_ANNOT:-${PROJECT_ROOT}/data/external/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl}"
CONDA_ENV="${SCENIC_CONDA_ENV:-pyscenic}"

N_WORKERS="${SLURM_CPUS_PER_TASK:-${SCENIC_WORKERS:-10}}"

mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}" "${TMP_BASE}"
export TMPDIR="${TMP_BASE}/${SLURM_JOB_ID:-manual}"
mkdir -p "${TMPDIR}"

# Activate the pySCENIC environment when conda is available. If conda is not
# available, the script assumes that python/pyscenic/Rscript are already on PATH.
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "${CONDA_ENV}"
elif [[ -n "${SCENIC_CONDA_SH:-}" && -f "${SCENIC_CONDA_SH}" ]]; then
  # shellcheck disable=SC1090
  source "${SCENIC_CONDA_SH}"
  conda activate "${CONDA_ENV}"
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

for cmd in Rscript python pyscenic; do
  command -v "${cmd}" >/dev/null 2>&1 || { echo "ERROR: ${cmd} not found on PATH" >&2; exit 1; }
done

for f in "${SEURAT_RDS}" "${TF_LIST}" "${RANKING_DB}" "${MOTIF_ANNOT}"; do
  [[ -s "${f}" ]] || { echo "ERROR: required input not found or empty: ${f}" >&2; exit 1; }
done

EXPR_CSV="${INPUT_DIR}/Stroma1_4_exprMat.csv"
LOOM_IN="${INPUT_DIR}/Stroma1_4_exprMat.loom"
ADJ="${OUTPUT_DIR}/adjacencies.tsv"
REG="${OUTPUT_DIR}/regulons.csv"
LOOM_OUT="${OUTPUT_DIR}/sce_SCENIC.loom"

echo "[INFO] Project root: ${PROJECT_ROOT}"
echo "[INFO] Workdir:      ${WORKDIR}"
echo "[INFO] Workers:      ${N_WORKERS}"
echo "[INFO] pySCENIC:     $(pyscenic --version 2>&1 || true)"

# Step 1. Extract raw Spatial counts from Stroma1_IGLC1 and Stroma4_IGHM.
Rscript --vanilla "${SCRIPT_DIR}/01_prepare_Stroma1_4_expr.R" \
  --seurat_rds "${SEURAT_RDS}" \
  --outdir "${INPUT_DIR}"

# Step 2. Convert the historical gene-by-spot CSV input to loom.
python "${SCRIPT_DIR}/02_csv_to_loom.py" \
  --expr "${EXPR_CSV}" \
  --out "${LOOM_IN}"

# Step 3. Infer transcription-factor/gene associations with GRNBoost2.
pyscenic grn \
  "${LOOM_IN}" \
  "${TF_LIST}" \
  --output "${ADJ}" \
  --num_workers "${N_WORKERS}" \
  --method grnboost2 \
  --seed 123

# Step 4. Motif enrichment / cisTarget pruning.
pyscenic ctx \
  "${ADJ}" \
  "${RANKING_DB}" \
  --annotations_fname "${MOTIF_ANNOT}" \
  --expression_mtx_fname "${LOOM_IN}" \
  --mode dask_multiprocessing \
  --output "${REG}" \
  --num_workers "${N_WORKERS}" \
  --mask_dropouts

# Step 5. AUCell regulon activity.
pyscenic aucell \
  "${LOOM_IN}" \
  "${REG}" \
  --output "${LOOM_OUT}" \
  --num_workers "${N_WORKERS}"

printf '%s\n' \
  "[DONE] pySCENIC core workflow finished." \
  "Next:" \
  "  Rscript --vanilla ${SCRIPT_DIR}/04_analyze_SCENIC_RSS_Stroma1_4.R --loom ${LOOM_OUT} --metadata ${INPUT_DIR}/Stroma1_4_metadata.csv --outdir ${OUTPUT_DIR}" \
  "  Rscript --vanilla ${SCRIPT_DIR}/05_plot_SCENIC_spatial_Stroma1_4.R --seurat_rds ${SEURAT_RDS} --auc_csv ${OUTPUT_DIR}/SCENIC_regulon_AUC_Stroma1_4.csv --outdir ${OUTPUT_DIR}/spatial"
