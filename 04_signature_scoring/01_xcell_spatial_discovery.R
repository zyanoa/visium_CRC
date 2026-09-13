# ==============================================================================
# Script: 01_xcell_spatial_discovery.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose: Reproduce the discovery-cohort xCell analysis used for Figure 3a
#          while avoiding destructive replacement of the Spatial assay.
#
# IMPORTANT:
#   The numerical workflow that generated the original Figure 3a is preserved:
#     Spatial raw counts -> xCellAnalysis -> CSV round-trip ->
#     CreateSeuratObject(counts = t(result)) -> NormalizeData -> ScaleData ->
#     mean scaled score by dMMR/pMMR spatial region.
#   The only structural change is that raw xCell scores used for spatial plotting
#   are stored in a separate "xCell" assay instead of overwriting Spatial@data.
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
})

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

DISCOVERY_ST_RDS <- file.path(PROCESSED_DATA_DIR, "st_obj.rds")
DISCOVERY_OUTDIR <- file.path(RESULTS_DIR, "04_signature_scoring", "xcell_discovery")
dir.create(DISCOVERY_OUTDIR, recursive = TRUE, showWarnings = FALSE)

get_assay_matrix <- function(object, assay, layer) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

# ------------------------------------------------------------------------------
# 1. Load discovery cohort
# ------------------------------------------------------------------------------

message("Loading discovery cohort: ", DISCOVERY_ST_RDS)
st_obj <- readRDS(DISCOVERY_ST_RDS)

required_meta <- c("type", "region")
missing_meta <- setdiff(required_meta, colnames(st_obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

# The historical script used PTdmmr$type_clusters. Preserve that grouping if the
# column is present; otherwise reconstruct the same dMMR/pMMR_region label.
if ("type_clusters" %in% colnames(st_obj@meta.data)) {
  st_obj$xcell_type_region <- as.character(st_obj$type_clusters)
} else {
  st_obj$xcell_type_region <- paste0(st_obj$type, "_", st_obj$region)
}

# ------------------------------------------------------------------------------
# 2. Run xCell on the same Spatial raw-count matrix used originally
# ------------------------------------------------------------------------------

message("Running xCellAnalysis on Spatial raw counts...")
expr_mat <- get_assay_matrix(st_obj, assay = "Spatial", layer = "counts")

xcell_result_raw <- xCell::xCellAnalysis(expr_mat)
xcell_result_raw <- as.data.frame(xcell_result_raw, check.names = FALSE)

# Preserve the historical intermediate file.
xcell_csv <- file.path(DISCOVERY_OUTDIR, "xCell_scores.csv")
write.csv(xcell_result_raw, xcell_csv, quote = FALSE)

# Historical code wrote and re-read the CSV with read.csv defaults. That converts
# labels such as "Plasma cells" and "B-cells" to "Plasma.cells" and "B.cells".
# Reproduce that behavior so downstream Figure 3a column names stay unchanged.
result <- read.csv(xcell_csv, row.names = 1, check.names = TRUE)
result <- as.data.frame(result)

spot_ids <- colnames(st_obj)

# Normalize orientation to: rows = spots, columns = xCell cell types.
if (all(spot_ids %in% rownames(result))) {
  result_by_spot <- result[spot_ids, , drop = FALSE]
} else if (all(spot_ids %in% colnames(result))) {
  result_by_spot <- as.data.frame(t(as.matrix(result[, spot_ids, drop = FALSE])))
} else {
  n_row_match <- sum(spot_ids %in% rownames(result))
  n_col_match <- sum(spot_ids %in% colnames(result))
  stop(
    "xCell output could not be aligned to the spatial object. ",
    "Spot matches in rows=", n_row_match, ", columns=", n_col_match,
    ", expected=", length(spot_ids)
  )
}

message("Aligned xCell score table: ", nrow(result_by_spot), " spots x ",
        ncol(result_by_spot), " cell types")

# ------------------------------------------------------------------------------
# 3. Reproduce the original normalization/scaling used for Figure 3a
# ------------------------------------------------------------------------------

# This is intentionally equivalent to the historical code:
# CreateSeuratObject(counts = t(result), meta.data = PTdmmr@meta.data)
xcell_obj <- CreateSeuratObject(
  counts = t(as.matrix(result_by_spot)),
  meta.data = st_obj@meta.data
)

xcell_obj <- NormalizeData(xcell_obj, assay = "RNA", verbose = FALSE)
xcell_obj <- ScaleData(xcell_obj, assay = "RNA", verbose = FALSE)
xcell_obj <- xcell_obj[, spot_ids]

scaled_data <- get_assay_matrix(xcell_obj, assay = "RNA", layer = "scale.data")
scaled_data <- as.data.frame(t(as.matrix(scaled_data)))
scaled_data$type_region <- st_obj$xcell_type_region[match(rownames(scaled_data), colnames(st_obj))]

if (anyNA(scaled_data$type_region)) {
  stop("Failed to align xCell scaled scores with spatial-region metadata.")
}

xcell_avg <- scaled_data %>%
  group_by(type_region) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

mat <- as.data.frame(xcell_avg[, -1, drop = FALSE])
rownames(mat) <- xcell_avg$type_region

# Original Figure 3a displayed these region pairs. Stroma6/7 were not part of the
# historical xCell panel and are therefore not added here, preserving the figure.
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

# New repository-facing name plus the historical filename used by the original analysis.
write.csv(mat, file.path(DISCOVERY_OUTDIR, "xcell_region_mean_scores.csv"), quote = FALSE)
write.csv(mat, file.path(DISCOVERY_OUTDIR, "Figure2G_XCELL.csv"), quote = FALSE)
saveRDS(xcell_obj, file.path(DISCOVERY_OUTDIR, "xcell_scaled_seurat.rds"))

# ------------------------------------------------------------------------------
# 4. Figure 3a B-cell lineage heatmap
# ------------------------------------------------------------------------------

b_features <- intersect(
  c("Plasma.cells", "B.cells", "naive.B.cells", "Memory.B.cells", "pro.B.cells"),
  colnames(mat)
)

if (length(b_features) > 0) {
  b_mat <- as.matrix(mat[, b_features, drop = FALSE])
  split_by_region <- gsub("^(dMMR_|pMMR_)", "", rownames(b_mat))
  split_by_region <- factor(split_by_region, levels = unique(split_by_region))

  pdf(file.path(DISCOVERY_OUTDIR, "xcell_Bcell_heatmap.pdf"), width = 8, height = 10)
  draw(Heatmap(
    b_mat,
    name = "xCell Score",
    col = colorRamp2(c(min(b_mat, na.rm = TRUE), max(b_mat, na.rm = TRUE)), c("white", "#CC0033")),
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
# 5. T-cell lineage heatmap used in the same historical workflow
# ------------------------------------------------------------------------------

t_features <- intersect(
  c(
    "CD4..memory.T.cells", "CD4..naive.T.cells", "CD4..T.cells", "CD4..Tcm",
    "CD8..naive.T.cells", "CD8..T.cells", "CD8..Tcm", "Tgd.cells",
    "Th1.cells", "Th2.cells", "Tregs"
  ),
  colnames(mat)
)

if (length(t_features) > 0) {
  t_mat <- as.matrix(mat[, t_features, drop = FALSE])
  split_by_region <- gsub("^(dMMR_|pMMR_)", "", rownames(t_mat))
  split_by_region <- factor(split_by_region, levels = unique(split_by_region))

  pdf(file.path(DISCOVERY_OUTDIR, "xcell_Tcell_heatmap.pdf"), width = 10, height = 10)
  draw(Heatmap(
    t_mat,
    name = "xCell Score",
    col = colorRamp2(c(min(t_mat, na.rm = TRUE), max(t_mat, na.rm = TRUE)), c("white", "#CC0033")),
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
# 6. Attach RAW xCell scores as a separate assay for spatial plotting
# ------------------------------------------------------------------------------

# Old GitHub code replaced Spatial@data with xCell scores. That did not change the
# xCell numbers, but it corrupted the meaning of the Spatial assay. Store the same
# raw xCell scores in their own assay instead.
xcell_feature_by_spot <- t(as.matrix(result_by_spot))
plot_obj <- st_obj
plot_obj[["xCell"]] <- CreateAssayObject(data = xcell_feature_by_spot)
DefaultAssay(plot_obj) <- "xCell"
saveRDS(plot_obj, file.path(DISCOVERY_OUTDIR, "st_obj_with_xCell_assay.rds"))

# ------------------------------------------------------------------------------
# 7. TLS-region Treg spatial plots
# ------------------------------------------------------------------------------

if ("Tregs" %in% rownames(plot_obj[["xCell"]])) {
  tls_obj <- subset(plot_obj, subset = region == "Stroma4_IGHM")
  DefaultAssay(tls_obj) <- "xCell"

  for (sid in unique(tls_obj$SectionID)) {
    p <- SpatialFeaturePlot(
      object = tls_obj,
      images = sid,
      features = "Tregs",
      pt.size.factor = 1.3,
      slot = "data"
    ) +
      scale_fill_gradientn(
        colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
        name = "xCell Score"
      ) +
      ggtitle(paste0(sid, " - Tregs")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(
      file.path(DISCOVERY_OUTDIR, paste0(sid, "_Tregs.pdf")),
      p, width = 6, height = 5
    )
  }
}

# ------------------------------------------------------------------------------
# 8. CR2 vs Tregs correlation in Stroma4_IGHM
# ------------------------------------------------------------------------------

if ("Tregs" %in% colnames(result_by_spot)) {
  tls_expr_obj <- subset(st_obj, subset = region == "Stroma4_IGHM")
  DefaultAssay(tls_expr_obj) <- "Spatial"
  tls_expr_obj <- NormalizeData(
    tls_expr_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  expression_CR2 <- FetchData(tls_expr_obj, vars = "CR2")
  treg_scores <- result_by_spot[colnames(tls_expr_obj), "Tregs"]

  cor_df <- data.frame(
    CR2 = expression_CR2$CR2,
    Tregs = as.numeric(treg_scores)
  )

  cor_test <- cor.test(cor_df$CR2, cor_df$Tregs, method = "spearman")

  cor_plot <- ggscatter(
    cor_df,
    x = "CR2",
    y = "Tregs",
    add = "reg.line",
    conf.int = TRUE,
    add.params = list(color = "#CC0033", fill = "lightgray")
  ) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_classic(base_size = 14) +
    labs(
      title = "TLS Region: CR2 Expression vs Tregs Score",
      subtitle = paste0(
        "Spearman rho = ", round(unname(cor_test$estimate), 3),
        ", p = ", format.pval(cor_test$p.value, digits = 3)
      ),
      x = "CR2 Expression (log-normalized)",
      y = "Tregs xCell Score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  ggsave(
    file.path(DISCOVERY_OUTDIR, "CR2_Tregs_correlation.pdf"),
    cor_plot, width = 6, height = 6
  )

  write.table(
    data.frame(
      n_spots = nrow(cor_df),
      spearman_rho = unname(cor_test$estimate),
      p_value = cor_test$p.value
    ),
    file.path(DISCOVERY_OUTDIR, "CR2_Tregs_correlation_stats.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(DISCOVERY_OUTDIR, "sessionInfo_xCell_spatial_discovery.txt")
)

message("xCell discovery-cohort analysis completed. Outputs: ", DISCOVERY_OUTDIR)
