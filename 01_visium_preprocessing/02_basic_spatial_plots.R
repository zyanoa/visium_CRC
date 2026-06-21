# ==============================================================================
# Script: 02_basic_spatial_plots.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   Generate basic QC and overview plots from the annotated Seurat object.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(RColorBrewer)
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

input_rds <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
marker_file <- file.path(RESULTS_DIR, "01_visium_preprocessing", "spatial_cluster_markers.csv")
outdir <- file.path(RESULTS_DIR, "01_visium_preprocessing", "basic_spatial_plots")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

st_obj <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# 2. Global plotting settings
# ------------------------------------------------------------------------------

section_levels <- c(
  "PT9_1", "PT9_2", "PT9_3",
  "PT18_1", "PT18_2", "PT18_3",
  "PT32_1", "PT32_2", "PT32_3",
  "PT34_1", "PT34_2",
  "PT55_1", "PT55_2", "PT55_3",
  "PT57_1", "PT57_2", "PT57_3", "PT57_4"
)

st_obj$SectionID <- factor(st_obj$SectionID, levels = section_levels)

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

Idents(st_obj) <- factor(Idents(st_obj), levels = cluster_levels)
st_obj$region <- factor(Idents(st_obj), levels = cluster_levels)

cluster_cols <- c(
  "Cancer1_CEACAM5" = "#e6194b",
  "Cancer2_KRT8" = "#911eb4",
  "Normal_PIGR" = "#fabebe",
  "Stroma1_IGLC1" = "#3cb44b",
  "Stroma2_MYL9+ACTG2" = "#4363d8",
  "Stroma3_APOE" = "#46f0f0",
  "Stroma4_IGHM" = "#ffe119",
  "Stroma5_CXCL8" = "#f58231",
  "Stroma6_MGP/MYH11" = "#f032e6",
  "Stroma7_FCGBP/COL3A1" = "#bcf60c"
)

section_cols <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4",
  "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff",
  "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1"
)
names(section_cols) <- section_levels

theme_pub <- theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# ------------------------------------------------------------------------------
# 3. QC violin plots
# ------------------------------------------------------------------------------

p_ncount <- VlnPlot(
  st_obj,
  features = "nCount_Spatial",
  group.by = "SectionID",
  pt.size = 0,
  log = TRUE
) +
  labs(x = NULL, y = "UMI counts (log scale)", title = "nCount_Spatial") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(outdir, "FigS_QC_nCount_Spatial.pdf"),
  plot = p_ncount,
  width = 14,
  height = 8
)

p_nfeature <- VlnPlot(
  st_obj,
  features = "nFeature_Spatial",
  group.by = "SectionID",
  pt.size = 0,
  log = TRUE
) +
  labs(x = NULL, y = "Detected genes (log scale)", title = "nFeature_Spatial") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(outdir, "FigS_QC_nFeature_Spatial.pdf"),
  plot = p_nfeature,
  width = 14,
  height = 8
)

# ------------------------------------------------------------------------------
# 4. UMAP overview plots
# ------------------------------------------------------------------------------

p_umap_section <- DimPlot(
  st_obj,
  reduction = "umap",
  group.by = "SectionID",
  cols = section_cols,
  raster = FALSE
) +
  labs(title = "UMAP by section") +
  theme_pub

ggsave(
  filename = file.path(outdir, "FigS_UMAP_by_SectionID.pdf"),
  plot = p_umap_section,
  width = 9,
  height = 6
)

p_umap_cluster <- DimPlot(
  st_obj,
  reduction = "umap",
  cols = cluster_cols,
  raster = FALSE
) +
  labs(title = "UMAP by spatial cluster") +
  theme_pub

ggsave(
  filename = file.path(outdir, "Fig2_UMAP_by_cluster.pdf"),
  plot = p_umap_cluster,
  width = 9,
  height = 6
)

# ------------------------------------------------------------------------------
# 5. Cluster composition by MMR status
# ------------------------------------------------------------------------------

spa_freq_df <- as.data.frame(
  prop.table(table(Idents(st_obj), st_obj$type), margin = 2)
)
colnames(spa_freq_df) <- c("Cluster", "MMR_status", "Fraction")

p_prop <- ggplot(spa_freq_df, aes(x = MMR_status, y = Fraction * 100, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.65, linewidth = 0.4) +
  scale_fill_manual(values = cluster_cols) +
  labs(x = NULL, y = "Percentage (%)", title = "Cluster composition by MMR status") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

ggsave(
  filename = file.path(outdir, "Fig2_cluster_composition_barplot.pdf"),
  plot = p_prop,
  width = 8,
  height = 6
)

# ------------------------------------------------------------------------------
# 6. Spatial cluster plots for each section
# ------------------------------------------------------------------------------

for (sample_id in levels(st_obj$SectionID)) {
  message("Generating spatial cluster plot for: ", sample_id)

  p_spatial <- SpatialDimPlot(
    object = st_obj,
    images = sample_id,
    cols = cluster_cols,
    pt.size.factor = 1.15,
    label = FALSE
  ) +
    labs(title = sample_id) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )

  ggsave(
    filename = file.path(outdir, paste0("Spatial_cluster_", sample_id, ".pdf")),
    plot = p_spatial,
    width = 6,
    height = 5
  )
}

# ------------------------------------------------------------------------------
# 7. Marker dot plot
# ------------------------------------------------------------------------------

marker_features <- c(
  "CEACAM5", "KRT8", "PIGR", "IGLC1",
  "MYL9", "ACTG2", "APOE", "IGHM",
  "CXCL8", "MGP", "MYH11", "FCGBP", "COL3A1"
)

p_dot <- DotPlot(st_obj, features = marker_features) +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradientn(
    colours = rev(c(
      "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
      "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
    ))
  ) +
  labs(title = "Canonical marker genes", x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  filename = file.path(outdir, "FigS_marker_dotplot.pdf"),
  plot = p_dot,
  width = 9,
  height = 6
)

# ------------------------------------------------------------------------------
# 8. Marker heatmap
# ------------------------------------------------------------------------------

diff.exp <- read.csv(marker_file)

top8 <- diff.exp %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 8) %>%
  ungroup()

p_heatmap <- DoHeatmap(
  object = subset(st_obj, downsample = 100),
  features = unique(top8$gene),
  assay = "SCT",
  group.colors = cluster_cols
) + NoLegend()

ggsave(
  filename = file.path(outdir, "Fig2_marker_heatmap.pdf"),
  plot = p_heatmap,
  width = 18,
  height = 12
)

# ------------------------------------------------------------------------------
# 9. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_02_basic_spatial_plots.txt")
)

message("All basic spatial plots have been generated successfully.")