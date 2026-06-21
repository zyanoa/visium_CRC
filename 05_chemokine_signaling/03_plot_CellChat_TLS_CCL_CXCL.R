#!/usr/bin/env Rscript

# Project: CRC_MMR_spatial_TLS
# Script: 03_plot_CellChat_TLS_CCL_CXCL.R
# Purpose:
#   Generate CellChat visualizations for TLS-centered chemokine signaling.
#   This script supports Figure 4F-G:
#     - Figure 4F-like rankNet panels for Stroma4->Stroma4, Stroma1->Stroma4,
#       and Stroma3->Stroma4 signaling.
#     - Figure 4G-like bubble plot for selected CCL/CXCL ligand-receptor pairs.

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    merged_cellchat = Sys.getenv("CRC_MMR_ST_CELLCHAT", unset = "results/05_chemokine_signaling/CellChat/cellchat_MMR_merged.rds"),
    output_dir = Sys.getenv("CRC_MMR_ST_RESULTS", unset = "results/05_chemokine_signaling/CellChat/Figure4FG")
  )
  for (a in args) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2]
  }
  opt
}

opt <- parse_args()
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

cellchat <- readRDS(opt$merged_cellchat)

# Numeric cluster indices follow the manuscript cluster order:
# 4 = Stroma1_IGLC1; 6 = Stroma3_APOE; 7 = Stroma4_IGHM / TLS.
ranknet_pairs <- list(
  Stroma4_to_Stroma4 = list(source = 7, target = 7),
  Stroma1_to_Stroma4 = list(source = 4, target = 7),
  Stroma3_to_Stroma4 = list(source = 6, target = 7)
)

for (nm in names(ranknet_pairs)) {
  source_id <- ranknet_pairs[[nm]]$source
  target_id <- ranknet_pairs[[nm]]$target

  p_stacked <- rankNet(
    cellchat,
    mode = "comparison",
    measure = "weight",
    sources.use = source_id,
    targets.use = target_id,
    stacked = TRUE,
    do.stat = TRUE,
    font.size = 12
  ) + ggtitle(nm)

  p_grouped <- rankNet(
    cellchat,
    mode = "comparison",
    measure = "weight",
    sources.use = source_id,
    targets.use = target_id,
    stacked = FALSE,
    do.stat = TRUE,
    font.size = 12
  ) + ggtitle(paste0(nm, " grouped"))

  ggsave(
    filename = file.path(opt$output_dir, paste0(nm, "_rankNet_stacked.pdf")),
    plot = p_stacked,
    width = 5,
    height = 6
  )

  ggsave(
    filename = file.path(opt$output_dir, paste0(nm, "_rankNet_grouped.pdf")),
    plot = p_grouped,
    width = 5,
    height = 6
  )

  ggsave(
    filename = file.path(opt$output_dir, paste0(nm, "_rankNet_combined.pdf")),
    plot = p_stacked + p_grouped,
    width = 10,
    height = 6
  )
}

# Selected CCL/CXCL ligand-receptor pairs used for the manuscript bubble plot.
pairLR <- c(
  "CCL19 - CCR7",
  "CCL20 - CCR6",
  "CCL21 - CCR7",
  "CCL3 - CCR5",
  "CCL4 - CCR5",
  "CCL5 - CCR5",
  "CXCL10 - CXCR3",
  "CXCL11 - CXCR3",
  "CXCL12 - CXCR4",
  "CXCL13 - CXCR3",
  "CXCL13 - CXCR5",
  "CXCL16 - CXCR6",
  "CXCL9 - CXCR3"
)
pairLR <- data.frame(interaction_name = gsub(" - ", "_", pairLR))

p_bubble <- netVisual_bubble(
  cellchat,
  sources.use = c(4, 6, 7),
  targets.use = 7,
  comparison = c(1, 2),
  angle.x = 45,
  pairLR.use = pairLR
)

ggsave(
  filename = file.path(opt$output_dir, "TLS_CCL_CXCL_selected_LR_bubble.pdf"),
  plot = p_bubble,
  width = 7,
  height = 5.5
)

message("Done. CellChat Figure 4F-G plots written to: ", opt$output_dir)
