# Spatial transcriptomic analysis of dMMR and pMMR colorectal cancer

This repository contains custom analysis scripts accompanying a spatial transcriptomics study of mismatch repair-deficient (dMMR) and mismatch repair-proficient (pMMR) colorectal cancer. The study investigates how spatially organized immune ecosystems, tertiary lymphoid structure maturation, plasma-cell-enriched niches, and chemokine-associated cellular interactions distinguish dMMR from pMMR tumors and relate to immunotherapy response and patient outcome.

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

## Data availability

Exceptions: Raw data can't be directly provided at the moment. Processed data will be available in proper data repositories after publication.
