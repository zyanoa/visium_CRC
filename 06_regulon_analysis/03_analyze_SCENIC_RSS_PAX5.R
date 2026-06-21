#!/usr/bin/env Rscript
# Analyze pySCENIC output and highlight PAX5 regulon specificity in Stroma3/Stroma4 groups.

suppressPackageStartupMessages({
  library(optparse)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(SummarizedExperiment)
})

option_list <- list(
  make_option("--loom", type = "character", default = "results/06_regulon_analysis/SCENIC_Stroma3_4/sce_SCENIC.loom"),
  make_option("--seurat_rds", type = "character", help = "Seurat object containing spatial metadata."),
  make_option("--output_dir", type = "character", default = "results/06_regulon_analysis"),
  make_option("--type_col", type = "character", default = "type", help = "Column containing dMMR/pMMR labels."),
  make_option("--region_col", type = "character", default = "region", help = "Column containing Stroma3_APOE/Stroma4_IGHM labels."),
  make_option("--target_regulon", type = "character", default = "PAX5(+)")
)
args <- parse_args(OptionParser(option_list = option_list))

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

sce_SCENIC <- open_loom(args$loom)
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = "RegulonsAUC")

obj <- readRDS(args$seurat_rds)
obj <- subset(obj, subset = .data[[args$region_col]] %in% c("Stroma3_APOE", "Stroma4_IGHM"))

obj$type_region <- paste0(obj@meta.data[[args$type_col]], "_", obj@meta.data[[args$region_col]])
target_groups <- c(
  "dMMR_Stroma3_APOE",
  "pMMR_Stroma3_APOE",
  "dMMR_Stroma4_IGHM",
  "pMMR_Stroma4_IGHM"
)
obj$type_region <- factor(obj$type_region, levels = target_groups)

common_cells <- intersect(colnames(obj), colnames(regulonAUC))
if (length(common_cells) == 0) stop("No matched cells/spots between Seurat object and SCENIC AUC matrix.")
obj <- obj[, common_cells]
sub_regulonAUC <- regulonAUC[, common_cells]

cellTypes <- data.frame(celltype = obj$type_region, row.names = colnames(obj))
rss <- calcRSS(
  AUC = getAUC(sub_regulonAUC),
  cellAnnotation = cellTypes[colnames(sub_regulonAUC), "celltype"]
)
rss <- na.omit(rss) |> as.data.frame()
write.csv(rss, file.path(args$output_dir, "SCENIC_RSS_Stroma3_4.csv"))

rss_df <- rss
if (!all(target_groups %in% colnames(rss_df))) {
  rss_df <- as.data.frame(t(as.matrix(rss_df)))
}
if (!all(target_groups %in% colnames(rss_df))) {
  stop("Target groups were not found in RSS matrix. Check metadata labels and row/column orientation.")
}

rss_df <- rss_df[, target_groups, drop = FALSE]
rss_df$Regulon <- rownames(rss_df)

rss_long <- rss_df |>
  pivot_longer(cols = all_of(target_groups), names_to = "Group", values_to = "RSS") |>
  group_by(Group) |>
  arrange(desc(RSS), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup()

top_df <- rss_long |>
  group_by(Group) |>
  slice_max(order_by = RSS, n = 3, with_ties = FALSE) |>
  ungroup() |>
  filter(Regulon != args$target_regulon)

highlight_df <- rss_long |> filter(Regulon == args$target_regulon)
if (nrow(highlight_df) == 0) {
  stop(paste0("Regulon not found: ", args$target_regulon,
              "\nUse grep('PAX5', rownames(rss), value = TRUE) to inspect available names."))
}

rss_long$Group <- factor(rss_long$Group, levels = target_groups)
top_df$Group <- factor(top_df$Group, levels = target_groups)
highlight_df$Group <- factor(highlight_df$Group, levels = target_groups)

tf_cols <- setNames("#D73027", args$target_regulon)

p <- ggplot(rss_long, aes(x = rank, y = RSS)) +
  geom_line(color = "grey65", linewidth = 0.8) +
  geom_point(data = top_df, color = "grey35", size = 1.8) +
  geom_text_repel(
    data = top_df,
    aes(label = Regulon),
    color = "grey25",
    size = 3.5,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = "grey60",
    max.overlaps = Inf
  ) +
  geom_point(data = highlight_df, aes(color = Regulon), size = 3.2) +
  geom_text_repel(
    data = highlight_df,
    aes(label = Regulon, color = Regulon),
    size = 4.8,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "black",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_wrap(~Group, scales = "free_y", nrow = 1) +
  scale_color_manual(values = tf_cols) +
  labs(x = "rank", y = "Regulon specificity score", title = "RSS ranking") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(args$output_dir, "PAX5_highlight_RSS_rank_4groups.pdf"), p, width = 6, height = 5)
ggsave(file.path(args$output_dir, "PAX5_highlight_RSS_rank_4groups.png"), p, width = 11, height = 6, dpi = 300)

# Also save AUC matrix for downstream spatial visualization.
auc_mat <- t(assay(sub_regulonAUC))
write.csv(auc_mat, file.path(args$output_dir, "SCENIC_regulon_AUC_matrix.csv"))
