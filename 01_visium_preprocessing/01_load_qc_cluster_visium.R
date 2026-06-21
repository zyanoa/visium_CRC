# ==============================================================================
# Script: 01_visium_qc_clustering.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Load 10x Visium spatial transcriptomics data
#   2. Perform QC metric calculation
#   3. Normalize data with SCTransform
#   4. Integrate sections and remove batch effects with Harmony
#   5. Perform clustering and marker detection
#   6. Annotate spatial clusters
#   7. Save the final Seurat object for downstream analyses
#
# Author: Zhengyang Zhao
# Date: 2026-06-21
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(RColorBrewer)
  library(harmony)
  library(ggpubr)
  library(pheatmap)
  library(grid)
  library(scRNAtoolVis)
  library(future)
})

rm(list = ls())
options(future.globals.maxSize = 20 * 1024^3)

# ------------------------------------------------------------------------------
# 1. Global settings
# ------------------------------------------------------------------------------

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

outdir <- file.path(RESULTS_DIR, "01_visium_preprocessing")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

SectionID <- c(
  paste0("PT18_", 1:3),
  paste0("PT32_", 1:3),
  paste0("PT34_", 1:2),
  paste0("PT57_", 1:4),
  paste0("PT55_", 1:3),
  paste0("PT9_", 1:3)
)

DataDir <- c(
  file.path(RAW_DATA_DIR, "PT18", "5CloupeFile", SectionID[1:3]),
  file.path(RAW_DATA_DIR, "PT32", "5CloupeFile", SectionID[4:6]),
  file.path(RAW_DATA_DIR, "PT34", "5CloupeFile", SectionID[7:8]),
  file.path(RAW_DATA_DIR, "PT57", "5CloupeFile", SectionID[9:12]),
  file.path(RAW_DATA_DIR, "PT55", "5CloupeFile", SectionID[13:15]),
  file.path(RAW_DATA_DIR, "PT9", "5CloupeFile", SectionID[16:18])
)

dmmr_ids <- c("PT18", "PT32", "PT9")
pmmr_ids <- c("PT34", "PT55", "PT57")

cluster_annotation_map <- c(
  "0" = "Cancer1_CEACAM5",
  "1" = "Stroma1_IGLC1",
  "2" = "Stroma2_MYL9+ACTG2",
  "3" = "Normal_PIGR",
  "4" = "Stroma3_APOE",
  "5" = "Stroma4_IGHM",
  "6" = "Stroma5_CXCL8",
  "7" = "Cancer2_KRT8",
  "8" = "Stroma6_MGP/MYH11",
  "9" = "Stroma7_FCGBP/COL3A1"
)

cluster_levels <- c(
  "Cancer1_CEACAM5",
  "Cancer2_KRT8",
  "Normal_PIGR",
  "Stroma1_IGLC1",
  "Stroma2_MYL9+ACTG2",
  "Stroma3_APOE",
  "Stroma4_IGHM",
  "Stroma5_CXCL8",
  "Stroma6_MGP/MYH11",
  "Stroma7_FCGBP/COL3A1"
)

# ------------------------------------------------------------------------------
# 2. Load Visium data and create Seurat objects
# ------------------------------------------------------------------------------

stolist <- vector("list", length = length(SectionID))
names(stolist) <- SectionID

for (i in seq_along(SectionID)) {
  message("Loading section: ", SectionID[i])

  tmp_obj <- Load10X_Spatial(
    data.dir = DataDir[i],
    slice = SectionID[i],
    filename = "filtered_feature_bc_matrix.h5"
  )

  tmp_obj <- RenameCells(tmp_obj, add.cell.id = SectionID[i])
  tmp_obj$SectionID <- SectionID[i]

  tmp_obj[["percent.mt"]] <- PercentageFeatureSet(tmp_obj, pattern = "^MT-")
  tmp_obj[["percent.rb"]] <- PercentageFeatureSet(tmp_obj, pattern = "^RP[LS]")

  stolist[[SectionID[i]]] <- tmp_obj
}

# ------------------------------------------------------------------------------
# 3. SCTransform normalization for each section
# ------------------------------------------------------------------------------

for (id in SectionID) {
  message("SCTransform: ", id)
  stolist[[id]] <- SCTransform(
    stolist[[id]],
    assay = "Spatial",
    verbose = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Merge all sections into one Seurat object
# ------------------------------------------------------------------------------

st_obj_raw <- merge(
  x = stolist[[1]],
  y = stolist[2:length(stolist)],
  add.cell.ids = names(stolist),
  project = "CRC_MMR_ST"
)

st_obj_raw@meta.data <- st_obj_raw@meta.data %>%
  mutate(
    patient_id = stringr::str_extract(SectionID, "PT\\d+"),
    type = case_when(
      patient_id %in% dmmr_ids ~ "dMMR",
      patient_id %in% pmmr_ids ~ "pMMR",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# 5. Integration features and dimensional reduction
# ------------------------------------------------------------------------------

VariableFeatures(st_obj_raw) <- unique(unlist(lapply(stolist, VariableFeatures)))

future::plan("multisession", workers = 8)

st_obj_raw <- RunPCA(
  st_obj_raw,
  assay = "SCT",
  verbose = FALSE
)

st_obj_raw <- RunHarmony(
  object = st_obj_raw,
  group.by.vars = "SectionID",
  assay.use = "SCT",
  reduction = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

st_obj_raw <- RunUMAP(
  object = st_obj_raw,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap"
)

st_obj_raw <- FindNeighbors(
  object = st_obj_raw,
  reduction = "harmony",
  dims = 1:20
)

st_obj_raw <- FindClusters(
  object = st_obj_raw,
  resolution = 0.1,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 6. Marker detection
# ------------------------------------------------------------------------------

st_obj_raw <- PrepSCTFindMarkers(
  object = st_obj_raw,
  assay = "SCT",
  verbose = TRUE
)

Idents(st_obj_raw) <- "seurat_clusters"

diff.exp <- FindAllMarkers(
  object = st_obj_raw,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.2
)

write.csv(
  diff.exp,
  file = file.path(outdir, "spatial_cluster_markers.csv"),
  row.names = FALSE
)

top_markers <- diff.exp %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

write.table(
  top_markers,
  file = file.path(outdir, "top5_markers_per_cluster.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Manual cluster annotation
# ------------------------------------------------------------------------------

st_obj <- st_obj_raw
Idents(st_obj) <- "seurat_clusters"

st_obj <- RenameIdents(st_obj, !!!cluster_annotation_map)
st_obj$region <- Idents(st_obj)

Idents(st_obj) <- factor(Idents(st_obj), levels = cluster_levels)
st_obj$region <- factor(st_obj$region, levels = cluster_levels)

# ------------------------------------------------------------------------------
# 8. Save final object
# ------------------------------------------------------------------------------

saveRDS(
  st_obj,
  file = file.path(outdir, "st_obj.rds")
)

saveRDS(
  st_obj_raw,
  file = file.path(outdir, "st_obj_raw_unannotated.rds")
)

message("Analysis completed successfully.")
message("Saved files:")
message("  - ", file.path(outdir, "st_obj.rds"))
message("  - ", file.path(outdir, "st_obj_raw_unannotated.rds"))
message("  - ", file.path(outdir, "spatial_cluster_markers.csv"))
message("  - ", file.path(outdir, "top5_markers_per_cluster.txt"))