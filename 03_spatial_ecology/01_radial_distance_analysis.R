# ==============================================================================
# Script: 01_radial_distance_analysis.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Compute radial distance from the tumor boundary for each spatial spot
#   2. Quantify inferred cell-type abundance as a function of distance
#   3. Export per-spot radial distance tables and smoothed distance plots
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Semla)
  library(ggplot2)
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

input_rds <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
cell2location_dir <- file.path(RESULTS_DIR, "02_celltype_annotation", "cell2location")
outdir <- file.path(RESULTS_DIR, "03_spatial_ecology", "radial_distance")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

st_obj <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# 2. Load and merge cell2location abundance matrices
# ------------------------------------------------------------------------------

sample_dirs <- list.dirs(cell2location_dir, full.names = TRUE, recursive = FALSE)
abundance_list <- list()

for (sub_dir in sample_dirs) {
  sample_name <- basename(sub_dir)
  abundance_file <- file.path(sub_dir, paste0(sample_name, ".spatial.deconvolution.csv"))

  if (file.exists(abundance_file)) {
    tmp_df <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
    abundance_list[[sample_name]] <- tmp_df
  }
}

if (length(abundance_list) == 0) {
  stop("No cell2location output files found.")
}

common_columns <- Reduce(intersect, lapply(abundance_list, colnames))
merged_abundance <- do.call(rbind, lapply(abundance_list, `[`, common_columns))
rownames(merged_abundance) <- gsub("\\.", "_", rownames(merged_abundance))
merged_abundance <- merged_abundance[colnames(st_obj), common_columns, drop = FALSE]

# ------------------------------------------------------------------------------
# 3. Define tumor ROI and compute radial distance
# ------------------------------------------------------------------------------

Idents(st_obj) <- "region"

section_ids <- unique(st_obj$SectionID)

radial_list <- list()

for (sid in section_ids) {
  message("Processing section: ", sid)

  slide_obj <- subset(st_obj, subset = SectionID == sid)

  # Tumor ROI defined according to manuscript: Cancer1_CEACAM5
  slide_obj$tumor_roi <- ifelse(Idents(slide_obj) == "Cancer1_CEACAM5", "tumor", "non_tumor")

  # Compute radial distance from the boundary of tumor ROI
  slide_obj <- RadialDistance(
    object = slide_obj,
    image_use = sid,
    group_by = "tumor_roi",
    selected_groups = "tumor",
    distance_name = "radial_distance"
  )

  radial_df <- FetchData(slide_obj, vars = c("SectionID", "region", "type", "radial_distance"))
  radial_df$spot_id <- rownames(radial_df)

  # Add cell2location abundance
  radial_df <- cbind(
    radial_df,
    merged_abundance[radial_df$spot_id, , drop = FALSE]
  )

  radial_list[[sid]] <- radial_df

  write.csv(
    radial_df,
    file = file.path(outdir, paste0(sid, "_radial_distance_table.csv")),
    row.names = FALSE
  )
}

radial_data <- bind_rows(radial_list)

write.csv(
  radial_data,
  file = file.path(outdir, "all_sections_radial_distance_table.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Example smoothed abundance-vs-distance plots
# ------------------------------------------------------------------------------

celltypes_to_plot <- intersect(
  c("B", "Plasma", "TCD8", "TCD4", "Macro"),
  colnames(radial_data)
)

for (ct in celltypes_to_plot) {
  p <- ggplot(radial_data, aes(x = radial_distance, y = .data[[ct]], color = type, fill = type)) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
    theme_classic(base_size = 14) +
    labs(
      x = "Distance from tumor boundary",
      y = "Inferred abundance",
      title = paste0(ct, " abundance across radial distance")
    )

  ggsave(
    filename = file.path(outdir, paste0("radial_distance_", ct, ".pdf")),
    plot = p,
    width = 6,
    height = 5
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_01_radial_distance_analysis.txt")
)