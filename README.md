# visium_CRC

This repository contains custom analysis scripts for a spatial transcriptomics study of mismatch repair-deficient (dMMR) and mismatch repair-proficient (pMMR) colorectal cancer, with a focus on tertiary lymphoid structure-associated immune microenvironments.

Raw sequencing data, large processed objects, and full sample-level metadata are not included in this repository. Software versions, sample information, and data sources are described in the manuscript and supplementary materials.

## Repository structure

```text
visium_CRC/
├── 00_utils/
├── 01_visium_preprocessing/
├── 02_celltype_annotation/
├── 03_spatial_ecology/
├── 04_signature_scoring/
├── 05_chemokine_signaling/
├── 06_regulon_analysis/
├── 07_external_validation/
└── signatures/
```

## Analysis modules

| Directory                  | Description                                                                  |
| -------------------------- | ---------------------------------------------------------------------------- |
| `00_utils/`                | Shared helper functions and project configuration                            |
| `01_visium_preprocessing/` | Visium preprocessing, quality control, clustering, and spatial visualization |
| `02_celltype_annotation/`  | Cell-type annotation and spatial deconvolution analyses                      |
| `03_spatial_ecology/`      | Spatial organization, niche, and microenvironment interaction analyses       |
| `04_signature_scoring/`    | Immune and stromal signature scoring in the discovery cohort                 |
| `05_chemokine_signaling/`  | Spatial chemokine-signaling analyses                                         |
| `06_regulon_analysis/`     | Transcription factor regulon activity analyses                               |
| `07_external_validation/`  | External cohort validation analyses                                          |
| `signatures/`              | Curated gene signatures used by selected validation scripts                  |

## Notes

The scripts are organized by analysis workflow rather than by figure panel. Local file paths should be adjusted in the project configuration file before running the analyses.

This repository is intended to accompany the manuscript and provide the custom code used for the reported analyses.
