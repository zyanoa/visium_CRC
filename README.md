# CRC MMR spatial TLS analysis code

This repository contains analysis scripts for a spatial transcriptomics study of tertiary lymphoid structure (TLS)-associated immune ecosystems in dMMR and pMMR colorectal cancer.

The repository is organized as a lightweight article-associated code release. Raw sequencing data, processed large objects, and full sample-level metadata are not included here. Software versions are described in the manuscript Methods section.

## Repository structure

```text
CRC_MMR_spatial_TLS/
├── README.md
├── LICENSE
├── .gitignore
├── 00_utils/
├── 01_visium_preprocessing/
├── 02_celltype_annotation/
├── 03_spatial_ecology/
├── 04_signature_scoring/
├── 05_chemokine_signaling/
├── 06_regulon_analysis/
├── 07_external_validation/
└── signatures/
    ├── marker_list.Rdata
    └── plasma_subtype_signatures.csv
```

## Analysis modules

| Module | Purpose |
|---|---|
| `00_utils/` | Shared helper functions and project path configuration |
| `01_visium_preprocessing/` | Visium loading, QC, clustering, and basic spatial visualization |
| `02_celltype_annotation/` | MIA and cell2location-based cell-type annotation |
| `03_spatial_ecology/` | Radial distance analysis, spatial niche analysis, and MISTy summaries |
| `04_signature_scoring/` | Discovery-cohort xCell scoring |
| `05_chemokine_signaling/` | COMMOT and CellChat chemokine-signaling analyses |
| `06_regulon_analysis/` | SCENIC/PAX5 regulon analysis |
| `07_external_validation/` | GSE236581 response validation, GSE39582/GSE17536 survival validation, and GSE13294 MSI/MSS plasma-cell validation |
| `signatures/` | Small signature files required by external validation scripts |

## Configuration

Default paths are defined in:

```text
00_utils/project_config.R
```

The scripts are written to avoid hard-coded personal paths. Local users should either edit `00_utils/project_config.R` or set environment variables before running the analysis.

Example:

```bash
export CRC_MMR_PROJECT_ROOT=/path/to/CRC_MMR_spatial_TLS
export CRC_MMR_DATA_ROOT=/path/to/local/data
export CRC_MMR_OUTPUT_ROOT=/path/to/local/output
```

## Typical running order

```bash
Rscript 01_visium_preprocessing/01_load_qc_cluster_visium.R
Rscript 01_visium_preprocessing/02_basic_spatial_plots.R

Rscript 02_celltype_annotation/01_prepare_MIA_markers.R
Rscript 02_celltype_annotation/02_run_MIA.R
sbatch  02_celltype_annotation/04_submit_cell2location_slurm.sh

Rscript 03_spatial_ecology/01_radial_distance_analysis.R
Rscript 03_spatial_ecology/02_spatial_niche_analysis.R
Rscript 03_spatial_ecology/04_run_misty_and_summarize.R
Rscript 03_spatial_ecology/05_misty_group_visualization.R

Rscript 04_signature_scoring/01_xcell_spatial_discovery.R
```

Representative chemokine-signaling analyses:

```bash
python 05_chemokine_signaling/01_run_COMMOT_CCL_CXCL_representative_samples.py \
  --input_dir data/processed/commot_h5ad \
  --output_dir results/05_chemokine_signaling/commot

Rscript 05_chemokine_signaling/02_run_CellChat_spatial_CCL_CXCL.R \
  --input_rds data/processed/st_object.rds \
  --output_dir results/05_chemokine_signaling/cellchat
```

SCENIC/PAX5 regulon analysis:

```bash
python 06_regulon_analysis/01_prepare_SCENIC_loom.py \
  --input_csv data/processed/scenic/Stroma3_4.csv \
  --output_loom results/06_regulon_analysis/sce.loom

sbatch 06_regulon_analysis/02_run_pySCENIC_slurm.sh

Rscript 06_regulon_analysis/03_analyze_SCENIC_RSS_PAX5.R \
  --loom results/06_regulon_analysis/sce_SCENIC.loom \
  --seurat_rds data/processed/st_object.rds \
  --output_dir results/06_regulon_analysis
```

External validation:

```bash
Rscript 07_external_validation/04_validate_GSE13294_MSI_plasma_programs.R \
  --expr_csv=data/external/GSE13294/GSE13294_exprSet_symbol.csv \
  --metadata_csv=data/external/GSE13294/GSE13294_metadata.csv \
  --plasma_signature_csv=signatures/plasma_subtype_signatures.csv \
  --run_xcell

Rscript 07_external_validation/01_validate_GSE236581_response.R \
  --seurat_rds=data/external/GSE236581/GSE236581_NT.rds \
  --marker_list=signatures/marker_list.Rdata

Rscript 07_external_validation/02_survival_KM_GSE39582_GSE17536.R \
  --gse39582_dir=data/external/GSE39582 \
  --marker_list=signatures/marker_list.Rdata

Rscript 07_external_validation/03_cox_forest_GSE17536.R
```

## Signature files

- `signatures/marker_list.Rdata` contains the final curated gene sets used for Figure 5 external validation. If the original object contains a generic `TLS` entry, the cleaned Figure 5 scripts exclude it. The internal `Stroma4_1` signature is renamed to `dMMR_TLS_program` in downstream Figure 5 scripts.
- Figure 3D uses the built-in xCell `Plasma cells` score.
- Figure 3E uses custom `B0_Plasma_IgA` and `B1_Plasma_IgG` gene sets stored in `signatures/plasma_subtype_signatures.csv`, derived from B-cell/plasma-cell subcluster marker genes.

## Notes

- Raw Visium, scRNA-seq, and external cohort expression matrices are not included.
- Full sample-level metadata should be reported in the manuscript supplementary tables.
- The cell2location scripts use `detection_alpha = 20`, matching the final analysis.
