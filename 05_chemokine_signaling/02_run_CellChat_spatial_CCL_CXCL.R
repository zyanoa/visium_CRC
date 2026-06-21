#!/usr/bin/env Rscript

# Project: CRC_MMR_spatial_TLS
# Script: 02_run_CellChat_spatial_CCL_CXCL.R
# Purpose:
#   Run spatial CellChat analysis for dMMR and pMMR CRC Visium sections.
#   This script focuses on secreted signaling and is intended to support the
#   Figure 4F-G analysis of TLS-centered CCL/CXCL communication.
#
# Main outputs:
#   1. results/05_chemokine_signaling/CellChat/cellchat_dMMR.rds
#   2. results/05_chemokine_signaling/CellChat/cellchat_pMMR.rds
#   3. results/05_chemokine_signaling/CellChat/cellchat_MMR_merged.rds
#
# Notes:
#   - Cluster labels must match the spatial cluster names used in the manuscript.
#   - The Stroma4_IGHM cluster is treated as the TLS-associated compartment.
#   - The script uses all dMMR or pMMR sections per group, not a single example section.
#   - The default CellChat settings follow the manuscript: spatial data, secreted
#     signaling, contact-dependent communication, and contact range = 100 um.

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(jsonlite)
  library(future)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    seurat_rds = Sys.getenv("CRC_MMR_ST_SEURAT_RDS", unset = "data/processed/st_obj.rds"),
    spaceranger_root = Sys.getenv("CRC_MMR_ST_RAW", unset = "data/raw/spaceranger"),
    output_dir = Sys.getenv("CRC_MMR_ST_RESULTS", unset = "results/05_chemokine_signaling/CellChat"),
    assay = "SCT",
    slot = "data",
    section_col = "SectionID",
    cluster_col = "seurat_clusters",
    spot_size_um = 55,
    interaction_range_um = 200,
    contact_range_um = 100,
    min_cells = 10,
    workers = 1
  )
  for (a in args) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2]
  }
  opt$spot_size_um <- as.numeric(opt$spot_size_um)
  opt$interaction_range_um <- as.numeric(opt$interaction_range_um)
  opt$contact_range_um <- as.numeric(opt$contact_range_um)
  opt$min_cells <- as.integer(opt$min_cells)
  opt$workers <- as.integer(opt$workers)
  opt
}

opt <- parse_args()
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Section grouping used in the study.
dmmr_sections <- c(
  "PT18_1", "PT18_2", "PT18_3",
  "PT32_1", "PT32_2", "PT32_3",
  "PT9_1",  "PT9_2",  "PT9_3"
)

pmmr_sections <- c(
  "PT34_1", "PT34_2",
  "PT55_1", "PT55_2", "PT55_3",
  "PT57_1", "PT57_2", "PT57_3", "PT57_4"
)

cluster_levels <- c(
  "Cancer1_CEACAM5",
  "Cancer2_KRT8",
  "Normal_PIGR",
  "Stroma1_IGLC1",
  "Stroma2_MYL9+ACTG2",
  "Stroma3_APOE",
  "Stroma4_IGHM",
  "Stroma5_CXCL8"
)

read_scalefactor <- function(section, spaceranger_root, spot_size_um) {
  patient_id <- strsplit(section, "_", fixed = TRUE)[[1]][1]
  candidate_paths <- c(
    file.path(spaceranger_root, patient_id, "5CloupeFile", section, "spatial", "scalefactors_json.json"),
    file.path(spaceranger_root, section, "spatial", "scalefactors_json.json")
  )
  scalefactors_path <- candidate_paths[file.exists(candidate_paths)][1]
  if (is.na(scalefactors_path)) {
    stop("Cannot find scalefactors_json.json for section: ", section,
         "\nTried:\n", paste(candidate_paths, collapse = "\n"))
  }
  sf <- jsonlite::fromJSON(scalefactors_path)
  conversion_factor <- spot_size_um / sf$spot_diameter_fullres
  data.frame(ratio = conversion_factor, tol = spot_size_um / 2)
}

build_cellchat_group <- function(seurat_obj, sections, group_name, opt) {
  message("Building CellChat object for ", group_name, " using sections: ", paste(sections, collapse = ", "))

  data_list <- list()
  meta_list <- list()
  locs_list <- list()
  factors_list <- list()

  for (section in sections) {
    object <- subset(seurat_obj, subset = .data[[opt$section_col]] == section)
    if (ncol(object) == 0) stop("No spots found for section: ", section)

    if (opt$cluster_col %in% colnames(object@meta.data)) {
      labels <- object@meta.data[[opt$cluster_col]]
    } else {
      labels <- Idents(object)
    }
    labels <- factor(as.character(labels), levels = cluster_levels)

    data_list[[section]] <- Seurat::GetAssayData(object, assay = opt$assay, slot = opt$slot)
    meta_list[[section]] <- data.frame(labels = labels, samples = section, group = group_name)

    locs <- Seurat::GetTissueCoordinates(
      object,
      scale = NULL,
      cols = c("imagerow", "imagecol"),
      image = section
    )
    locs_list[[section]] <- locs
    factors_list[[section]] <- read_scalefactor(section, opt$spaceranger_root, opt$spot_size_um)
  }

  data_input <- do.call(cbind, data_list)
  meta <- do.call(rbind, meta_list)
  rownames(meta) <- colnames(data_input)

  spatial_locs <- do.call(rbind, locs_list)
  rownames(spatial_locs) <- colnames(data_input)

  spatial_factors <- do.call(rbind, factors_list)

  meta$labels <- factor(meta$labels, levels = cluster_levels)

  cellchat <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial_locs,
    spatial.factors = spatial_factors
  )
  cellchat@idents <- factor(cellchat@meta$labels, levels = cluster_levels)

  CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # Manuscript-oriented setting: spatial communication with local contact-dependent constraints.
  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = FALSE,
    interaction.range = opt$interaction_range_um,
    scale.distance = NULL,
    contact.dependent = TRUE,
    contact.range = opt$contact_range_um
  )

  cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}

message("Loading Seurat object: ", opt$seurat_rds)
seurat_obj <- readRDS(opt$seurat_rds)

if (opt$workers > 1) {
  future::plan("multisession", workers = opt$workers)
} else {
  future::plan("sequential")
}

cellchat_dmmr <- build_cellchat_group(seurat_obj, dmmr_sections, "dMMR", opt)
cellchat_pmmr <- build_cellchat_group(seurat_obj, pmmr_sections, "pMMR", opt)

saveRDS(cellchat_dmmr, file.path(opt$output_dir, "cellchat_dMMR.rds"))
saveRDS(cellchat_pmmr, file.path(opt$output_dir, "cellchat_pMMR.rds"))

cellchat_merged <- mergeCellChat(
  list(dMMR = cellchat_dmmr, pMMR = cellchat_pmmr),
  add.names = c("dMMR", "pMMR")
)
saveRDS(cellchat_merged, file.path(opt$output_dir, "cellchat_MMR_merged.rds"))

message("Done. Outputs written to: ", opt$output_dir)
