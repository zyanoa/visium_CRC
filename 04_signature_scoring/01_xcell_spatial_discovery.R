# ==============================================================================
# Script: 01_xcell_spatial_discovery.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose: xCell analysis for the discovery spatial transcriptomics cohort
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(xCell)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
  library(rstatix)
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Discovery cohort (Spatial Transcriptomics)
DISCOVERY_ST_RDS <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
DISCOVERY_OUTDIR <- file.path(RESULTS_DIR, "04_signature_scoring", "xcell_discovery")

# Create output directory
dir.create(DISCOVERY_OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PART 1: DISCOVERY COHORT - SPATIAL TRANSCRIPTOMICS
# ==============================================================================

message("==========================================")
message("PART 1: Discovery Cohort Analysis")
message("==========================================")

# ------------------------------------------------------------------------------
# 1.1 Load discovery cohort data
# ------------------------------------------------------------------------------

message("Loading discovery cohort spatial transcriptomics data...")
st_obj <- readRDS(DISCOVERY_ST_RDS)

# ------------------------------------------------------------------------------
# 1.2 Run xCell on spatial count matrix
# ------------------------------------------------------------------------------

message("Running xCell analysis on spatial transcriptomics data...")
expr_mat <- st_obj@assays$Spatial@counts
xcell_scores <- xCell::xCellAnalysis(expr_mat)

write.csv(
  xcell_scores,
  file = file.path(DISCOVERY_OUTDIR, "xCell_scores.csv"),
  quote = FALSE
)

message("xCell scores saved to: ", DISCOVERY_OUTDIR)

# ------------------------------------------------------------------------------
# 1.3 Create Seurat object with xCell scores
# ------------------------------------------------------------------------------

message("Creating Seurat object with xCell scores...")
xcell_obj <- CreateSeuratObject(
  counts = t(xcell_scores),
  meta.data = st_obj@meta.data
)

xcell_obj <- NormalizeData(xcell_obj, assay = "RNA")
xcell_obj <- ScaleData(xcell_obj)
xcell_obj <- xcell_obj[, colnames(st_obj)]

saveRDS(
  xcell_obj,
  file = file.path(DISCOVERY_OUTDIR, "xcell_obj.rds")
)

# ------------------------------------------------------------------------------
# 1.4 Region-level summary statistics
# ------------------------------------------------------------------------------

message("Calculating region-level mean xCell scores...")
st_obj$type_region <- paste0(st_obj$type, "_", st_obj$region)

scaled_data <- xcell_obj[["RNA"]]@scale.data
scaled_data <- as.data.frame(t(scaled_data))
scaled_data$type_region <- st_obj$type_region

xcell_avg <- scaled_data %>%
  group_by(type_region) %>%
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE))

mat <- as.matrix(xcell_avg[, -1])
rownames(mat) <- xcell_avg$type_region

# Sort regions in biologically meaningful order
sorted_regions <- c(
  "dMMR_Cancer1_CEACAM5", "pMMR_Cancer1_CEACAM5",
  "dMMR_Cancer2_KRT8", "pMMR_Cancer2_KRT8",
  "dMMR_Normal_PIGR", "pMMR_Normal_PIGR",
  "dMMR_Stroma1_IGLC1", "pMMR_Stroma1_IGLC1",
  "dMMR_Stroma2_MYL9+ACTG2", "pMMR_Stroma2_MYL9+ACTG2",
  "dMMR_Stroma3_APOE", "pMMR_Stroma3_APOE",
  "dMMR_Stroma4_IGHM", "pMMR_Stroma4_IGHM",
  "dMMR_Stroma5_CXCL8", "pMMR_Stroma5_CXCL8"
)

valid_regions <- intersect(sorted_regions, rownames(mat))
mat <- mat[valid_regions, , drop = FALSE]

write.csv(
  mat,
  file = file.path(DISCOVERY_OUTDIR, "xcell_region_mean_scores.csv"),
  quote = FALSE
)

# ------------------------------------------------------------------------------
# 1.5 B-cell lineage heatmap
# ------------------------------------------------------------------------------

message("Generating B-cell lineage heatmap...")
b_features <- intersect(
  c("Plasma.cells", "B.cells", "naive.B.cells", "Memory.B.cells", "pro.B.cells"),
  colnames(mat)
)

if (length(b_features) > 0) {
  b_mat <- mat[, b_features, drop = FALSE]
  split_by_region <- gsub("^(dMMR_|pMMR_)", "", rownames(b_mat))
  split_by_region <- factor(split_by_region, levels = unique(split_by_region))

  pdf(file.path(DISCOVERY_OUTDIR, "xcell_Bcell_heatmap.pdf"), width = 8, height = 10)
  draw(Heatmap(
    as.matrix(b_mat),
    name = "xCell Score",
    col = colorRamp2(c(0, max(b_mat, na.rm = TRUE)), c("white", "#CC0033")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "black"),
    row_split = split_by_region,
    column_title = "B-cell Lineage Scores",
    row_title = "Region",
    heatmap_legend_param = list(title = "xCell Score")
  ))
  dev.off()
}

# ------------------------------------------------------------------------------
# 1.6 T-cell lineage heatmap
# ------------------------------------------------------------------------------

message("Generating T-cell lineage heatmap...")
t_features <- intersect(
  c(
    "CD4..memory.T.cells", "CD4..naive.T.cells", "CD4..T.cells", "CD4..Tcm",
    "CD8..naive.T.cells", "CD8..T.cells", "CD8..Tcm", "Tgd.cells",
    "Th1.cells", "Th2.cells", "Tregs"
  ),
  colnames(mat)
)

if (length(t_features) > 0) {
  t_mat <- mat[, t_features, drop = FALSE]
  split_by_region <- gsub("^(dMMR_|pMMR_)", "", rownames(t_mat))
  split_by_region <- factor(split_by_region, levels = unique(split_by_region))

  pdf(file.path(DISCOVERY_OUTDIR, "xcell_Tcell_heatmap.pdf"), width = 10, height = 10)
  draw(Heatmap(
    as.matrix(t_mat),
    name = "xCell Score",
    col = colorRamp2(c(0, max(t_mat, na.rm = TRUE)), c("white", "#CC0033")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "black"),
    row_split = split_by_region,
    column_title = "T-cell Lineage Scores",
    row_title = "Region",
    heatmap_legend_param = list(title = "xCell Score")
  ))
  dev.off()
}

# ------------------------------------------------------------------------------
# 1.7 Add xCell scores to ST object for spatial plotting
# ------------------------------------------------------------------------------

message("Adding xCell scores to ST object...")
plot_obj <- st_obj
plot_obj@assays$Spatial@data <- t(as.data.frame(xcell_scores))

# ------------------------------------------------------------------------------
# 1.8 TLS region spatial plots
# ------------------------------------------------------------------------------

message("Generating spatial plots for TLS region...")
tls_obj <- subset(plot_obj, subset = region == "Stroma4_IGHM")

features_to_plot <- c("Tregs")
features_to_plot <- intersect(features_to_plot, rownames(tls_obj))

for (feature_name in features_to_plot) {
  for (sid in unique(tls_obj$SectionID)) {
    p <- SpatialFeaturePlot(
      object = tls_obj,
      images = sid,
      features = feature_name,
      pt.size.factor = 1.3,
      slot = "data"
    ) +
      scale_fill_gradientn(
        colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
        name = "xCell Score"
      ) +
      ggtitle(paste0(sid, " - ", feature_name)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(
      filename = file.path(DISCOVERY_OUTDIR, paste0(sid, "_", feature_name, ".pdf")),
      plot = p,
      width = 6,
      height = 5
    )
  }
}

# ------------------------------------------------------------------------------
# 1.9 CR2 vs Tregs correlation in TLS spots
# ------------------------------------------------------------------------------

message("Analyzing CR2 vs Tregs correlation in TLS region...")
tls_expr_obj <- subset(st_obj, subset = region == "Stroma4_IGHM")
tls_expr_obj <- NormalizeData(tls_expr_obj, normalization.method = "LogNormalize", scale.factor = 10000)

if ("Tregs" %in% rownames(xcell_obj)) {
  expression_CR2 <- FetchData(tls_expr_obj, vars = "CR2")
  xcell_raw <- read.csv(
    file.path(DISCOVERY_OUTDIR, "xCell_scores.csv"),
    row.names = 1,
    check.names = FALSE
  )
  xcell_raw <- as.data.frame(t(xcell_raw))
  treg_scores <- xcell_raw[colnames(tls_expr_obj), "Tregs"]

  cor_df <- data.frame(
    CR2 = expression_CR2$CR2,
    Tregs = treg_scores
  )

  # Spearman correlation test
  cor_test <- cor.test(cor_df$CR2, cor_df$Tregs, method = "spearman")
  
  cor_plot <- ggscatter(
    cor_df,
    x = "CR2",
    y = "Tregs",
    add = "reg.line",
    conf.int = TRUE,
    add.params = list(color = "#CC0033", fill = "lightgray")
  ) +
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top"
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = "TLS Region: CR2 Expression vs Tregs Score",
      subtitle = paste0("Spearman ρ = ", round(cor_test$estimate, 3),
                        ", p = ", format.pval(cor_test$p.value, digits = 3)),
      x = "CR2 Expression (log-normalized)",
      y = "Tregs xCell Score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  ggsave(
    filename = file.path(DISCOVERY_OUTDIR, "CR2_Tregs_correlation.pdf"),
    plot = cor_plot,
    width = 6,
    height = 6
  )
  
  # Save correlation statistics
  sink(file.path(DISCOVERY_OUTDIR, "CR2_Tregs_correlation_stats.txt"))
  cat("CR2 vs Tregs Correlation Analysis (TLS Region)\n")
  cat("==============================================\n\n")
  cat("Number of spots:", nrow(cor_df), "\n")
  cat("Spearman correlation coefficient:", cor_test$estimate, "\n")
  cat("P-value:", cor_test$p.value, "\n")
  sink()
}


# ==============================================================================
# SESSION INFORMATION
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(DISCOVERY_OUTDIR, "sessionInfo_xCell_spatial_discovery.txt")
)

message("xCell discovery-cohort analysis completed. Outputs written to: ", DISCOVERY_OUTDIR)
