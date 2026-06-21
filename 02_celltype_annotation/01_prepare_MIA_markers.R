# ==============================================================================
# Script: 01_MIA_prepare_markers.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Load annotated spatial Seurat object and scRNA-seq reference object
#   2. Identify marker genes for scRNA-seq cell types
#   3. Identify marker genes for spatial regions
#   4. Export formatted marker tables for downstream MIA analysis
#
# Input:
#   - st_obj.rds
#   - scRNA reference object (Seurat)
#
# Output:
#   - sc_markers_full.csv
#   - spatial_markers_full.csv
#   - celltype_specific_for_MIA.csv
#   - region_specific_for_MIA.csv
#
# Author: Zhengyang Zhao
# Date: 2026-06-21
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()
options(future.globals.maxSize = 20 * 1024^3)

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

input_st_rds <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
input_sc_rds <- file.path(PROCESSED_DATA_DIR, "crc_scRNA_reference.rds")
outdir <- file.path(RESULTS_DIR, "02_celltype_annotation", "MIA")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

st_obj <- readRDS(input_st_rds)
cortex_sc <- readRDS(input_sc_rds)

# ------------------------------------------------------------------------------
# 2. Global settings
# ------------------------------------------------------------------------------

# Spatial region order used in MIA heatmap
region_levels <- c(
  "dMMR_Cancer1_CEACAM5",
  "pMMR_Cancer1_CEACAM5",
  "dMMR_Cancer2_KRT8",
  "pMMR_Cancer2_KRT8",
  "dMMR_Normal_PIGR",
  "pMMR_Normal_PIGR",
  "dMMR_Stroma1_IGLC1",
  "pMMR_Stroma1_IGLC1",
  "dMMR_Stroma2_MYL9+ACTG2",
  "pMMR_Stroma2_MYL9+ACTG2",
  "dMMR_Stroma3_APOE",
  "pMMR_Stroma3_APOE",
  "dMMR_Stroma4_IGHM",
  "pMMR_Stroma4_IGHM",
  "dMMR_Stroma5_CXCL8",
  "pMMR_Stroma5_CXCL8"
)

# Remove sparse or absent pMMR stromal clusters if needed
excluded_regions <- c(
  "pMMR_Stroma6_MGP/MYH11",
  "pMMR_Stroma7_FCGBP/COL3A1"
)

# scRNA marker thresholds
sc_logfc_cutoff <- 0.25
sc_min_pct <- 0.25
sc_delta_pct_cutoff <- 0.20
sc_padj_cutoff <- 0.05

# Spatial marker thresholds
sp_logfc_cutoff <- 0.20
sp_min_pct <- 0.10
sp_delta_pct_cutoff <- 0.10
sp_padj_cutoff <- 0.05

# ------------------------------------------------------------------------------
# 3. Prepare scRNA-seq markers
# ------------------------------------------------------------------------------

DefaultAssay(cortex_sc) <- "RNA"
Idents(cortex_sc) <- cortex_sc$clMidwayPr

sc.markers <- FindAllMarkers(
  object = cortex_sc,
  only.pos = TRUE,
  min.pct = sc_min_pct,
  logfc.threshold = sc_logfc_cutoff
)

sc.markers <- sc.markers %>%
  mutate(delta_pct = pct.1 - pct.2) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  sc.markers,
  file = file.path(outdir, "sc_markers_full.csv"),
  row.names = FALSE
)

sc.main.marker <- sc.markers %>%
  filter(
    avg_log2FC > sc_logfc_cutoff,
    p_val_adj < sc_padj_cutoff,
    delta_pct > sc_delta_pct_cutoff
  ) %>%
  arrange(cluster, desc(avg_log2FC))

celltype_specific <- sc.main.marker %>%
  dplyr::select(cluster, gene) %>%
  dplyr::rename(celltype = cluster)

write.csv(
  celltype_specific,
  file = file.path(outdir, "celltype_specific_for_MIA.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Prepare spatial region markers
# ------------------------------------------------------------------------------

st_obj$type_region <- paste0(st_obj$type, "_", Idents(st_obj))
Idents(st_obj) <- "type_region"

# Exclude very sparse regions if required by the MIA design
st_obj_mia <- subset(st_obj, idents = excluded_regions, invert = TRUE)

Idents(st_obj_mia) <- factor(Idents(st_obj_mia), levels = region_levels)

sp.markers <- FindAllMarkers(
  object = st_obj_mia,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = sp_min_pct,
  logfc.threshold = sp_logfc_cutoff
)

sp.markers <- sp.markers %>%
  mutate(delta_pct = pct.1 - pct.2) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(
  sp.markers,
  file = file.path(outdir, "spatial_markers_full.csv"),
  row.names = FALSE
)

region_main_marker <- sp.markers %>%
  filter(
    avg_log2FC > sp_logfc_cutoff,
    p_val_adj < sp_padj_cutoff,
    delta_pct > sp_delta_pct_cutoff
  ) %>%
  arrange(cluster, desc(avg_log2FC))

region_specific <- region_main_marker %>%
  dplyr::select(cluster, gene) %>%
  dplyr::rename(region = cluster)

write.csv(
  region_specific,
  file = file.path(outdir, "region_specific_for_MIA.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_01_MIA_prepare_markers.txt")
)

message("MIA marker preparation completed successfully.")
message("Saved files:")
message("  - ", file.path(outdir, "sc_markers_full.csv"))
message("  - ", file.path(outdir, "spatial_markers_full.csv"))
message("  - ", file.path(outdir, "celltype_specific_for_MIA.csv"))
message("  - ", file.path(outdir, "region_specific_for_MIA.csv"))