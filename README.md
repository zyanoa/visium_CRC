# visium_CRC

This repository contains custom analysis scripts for a spatial transcriptomics study of mismatch repair-deficient (dMMR) and mismatch repair-proficient (pMMR) colorectal cancer, with a focus on tertiary lymphoid structure (TLS)-associated immune ecosystems.

Raw sequencing data, large processed objects, and full sample-level metadata are not included in this repository. Software versions and sample information are described in the manuscript and supplementary tables.

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

| Directory                  | Description                                                            |
| -------------------------- | ---------------------------------------------------------------------- |
| `00_utils/`                | Shared helper functions and project configuration                      |
| `01_visium_preprocessing/` | Visium loading, quality control, clustering, and spatial visualization |
| `02_celltype_annotation/`  | MIA and cell2location-based cell-type annotation                       |
| `03_spatial_ecology/`      | Radial distance analysis, spatial niche analysis, and MISTy analysis   |
| `04_signature_scoring/`    | Discovery-cohort immune signature scoring                              |
| `05_chemokine_signaling/`  | COMMOT and CellChat chemokine-signaling analysis                       |
| `06_regulon_analysis/`     | SCENIC/PAX5 regulon analysis                                           |
| `07_external_validation/`  | External response, survival, and MSI/MSS validation analyses           |
| `signatures/`              | Small signature files required by external validation scripts          |

## Signature files

`signatures/marker_list.Rdata` contains curated gene sets used for external validation. In the cleaned scripts, the internal `Stroma4_1` signature is renamed to `dMMR_TLS_program`, and the generic `TLS` entry is excluded from the final Figure 5 analyses.

`signatures/plasma_subtype_signatures.csv` contains the custom `B0_Plasma_IgA` and `B1_Plasma_IgG` gene sets used for ssGSEA-based plasma-cell subtype validation. The Figure 3D plasma-cell score was estimated using the built-in xCell `Plasma cells` signature.

## Notes

The scripts are organized by analysis workflow rather than by figure panel. Local paths should be adjusted in `00_utils/project_config.R` before running the analyses.

The cell2location scripts use `detection_alpha = 20`, matching the final analysis.
