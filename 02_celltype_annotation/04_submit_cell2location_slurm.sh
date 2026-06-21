#!/bin/bash
#SBATCH --job-name=cell2location_crc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=16
#SBATCH --error=%j.err
#SBATCH --output=%j.out

set -euo pipefail

source ~/.bashrc
conda activate cell2location

# ------------------------------------------------------------------------------
# 1. Input settings
# ------------------------------------------------------------------------------

PROJECT_DIR="${CRC_MMR_ST_PROJECT:-$(pwd)}"
DATA_DIR="${CRC_MMR_ST_DATA:-${PROJECT_DIR}/data}"
RAW_BASE="${CRC_MMR_ST_RAW:-${DATA_DIR}/raw}"
RESULTS_DIR="${PROJECT_DIR}/results"

INPUT_SC="${DATA_DIR}/processed/sc_adata.h5ad"
SC_TYPE="10x"
PY_SCRIPT="${PROJECT_DIR}/02_celltype_annotation/03_run_cell2location.py"

OUTDIR_BASE="${RESULTS_DIR}/02_celltype_annotation/cell2location"

# ------------------------------------------------------------------------------
# 2. Patient-level sample definition
# ------------------------------------------------------------------------------

declare -A PATIENT_SECTION_COUNT=(
  [PT18]=3
  [PT32]=3
  [PT34]=2
  [PT55]=3
  [PT57]=4
  [PT9]=3
)

# ------------------------------------------------------------------------------
# 3. Run cell2location
# ------------------------------------------------------------------------------

for PATIENT_ID in PT18 PT32 PT34 PT55 PT57 PT9; do
    N_SECTIONS="${PATIENT_SECTION_COUNT[$PATIENT_ID]}"

    for IDX in $(seq 1 "${N_SECTIONS}"); do
        SAMPLE_NAME="${PATIENT_ID}_${IDX}"
        SP_INPUT="${RAW_BASE}/${PATIENT_ID}/5CloupeFile/${SAMPLE_NAME}"
        SAMPLE_OUTDIR="${OUTDIR_BASE}/${SAMPLE_NAME}"

        echo "[INFO] Running cell2location for ${SAMPLE_NAME}"
        echo "[INFO] Spatial input: ${SP_INPUT}"
        echo "[INFO] Output dir: ${SAMPLE_OUTDIR}"

        python "${PY_SCRIPT}" \
            --input_sc "${INPUT_SC}" \
            --sc_type "${SC_TYPE}" \
            --input_sp "${SP_INPUT}" \
            --outdir "${SAMPLE_OUTDIR}" \
            --sample "${SAMPLE_NAME}" \
            --detection_alpha 20
    done
done