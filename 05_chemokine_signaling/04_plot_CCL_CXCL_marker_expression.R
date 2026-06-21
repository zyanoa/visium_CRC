#!/usr/bin/env Rscript

# Project: CRC_MMR_spatial_TLS
# Script: 04_plot_CCL_CXCL_marker_expression.R
# Purpose:
#   Plot spatial expression and DotPlot summaries for key chemokine genes that
#   support the CellChat/COMMOT interpretation: CXCL13, CXCR5, CCL19, CCL21, CCR7.
#   This script supports the Figure 4H-I marker-expression panels.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    seurat_rds = Sys.getenv("CRC_MMR_ST_SEURAT_RDS", unset = "data/processed/st_obj.rds"),
    output_dir = Sys.getenv("CRC_MMR_ST_RESULTS", unset = "results/05_chemokine_signaling/CellChat/marker_expression"),
    sections = "PT9_3,PT55_2",
    assay = "SCT",
    slot = "data",
    max_cutoff = 2
  )
  for (a in args) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2]
  }
  opt$sections <- unlist(strsplit(opt$sections, ","))
  opt$max_cutoff <- as.numeric(opt$max_cutoff)
  opt
}

opt <- parse_args()
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$seurat_rds)
DefaultAssay(obj) <- opt$assay

features <- c("CXCL13", "CXCR5", "CCL19", "CCL21", "CCR7")

for (feature in features) {
  for (section in opt$sections) {
    p <- SpatialFeaturePlot(
      object = obj,
      images = section,
      features = feature,
      pt.size.factor = 1.3,
      alpha = c(0.6, 1),
      max.cutoff = opt$max_cutoff,
      slot = opt$slot,
      combine = TRUE
    ) +
      scale_fill_gradientn(
        colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
        limits = c(0, opt$max_cutoff)
      ) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = "right"
      ) +
      ggtitle(paste0(section, "_", feature))

    ggsave(
      filename = file.path(opt$output_dir, paste0(section, "_", feature, ".pdf")),
      plot = p,
      width = 6,
      height = 5
    )
  }
}

p_dot <- DotPlot(obj, features = c("CCL19", "CCR7", "CXCL13", "CXCR5")) +
  RotatedAxis() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(size = guide_legend("Percent Expression")) +
  scale_color_gradientn(
    colours = rev(c(
      "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
      "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
    ))
  )

ggsave(
  filename = file.path(opt$output_dir, "CCL_CXCL_marker_DotPlot.pdf"),
  plot = p_dot,
  width = 9,
  height = 6
)

message("Done. Marker-expression plots written to: ", opt$output_dir)
