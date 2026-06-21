# ==============================================================================
# Script: 02_niche_analysis.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Merge cell2location abundance matrices across sections
#   2. Perform ILR transformation of cell compositions
#   3. Identify recurrent spatial niches using SNN + Louvain clustering
#   4. Generate UMAP, niche summary, and niche abundance heatmaps
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(compositions)
  library(clustree)
  library(uwot)
  library(scran)
  library(cluster)
  library(RColorBrewer)
  library(dplyr)
  library(ggrastr)
  library(broom)
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
outdir <- file.path(RESULTS_DIR, "03_spatial_ecology", "spatial_niches")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

st_obj <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# 2. Colors and parameters
# ------------------------------------------------------------------------------

defined_cols <- c(
  "#B4DEF7", "#70AA43", "#FCD700", "#8150A1",
  "#FCA023", "#AFD242", "#595B5A", "#D1AF81",
  "#4B5CAA", "#EC2124", "#A86E3C", "#7F7782"
)

k_vect <- c(10, 20, 30)
selected_k <- 30
selected_resolution <- 0.5

celltype_order <- c(
  "Epi","TCD4", "TCD8","Tgd","NK","Mast","TZBTB16",
  "B", "Plasma","ILC","Peri","Endo","Fibro",
  "SmoothMuscle","Granulo","Mono","Macro","DC","Schwann"
)

# ------------------------------------------------------------------------------
# 3. Read and merge cell2location abundance matrices
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
  stop("No cell2location abundance files found.")
}

common_columns <- Reduce(intersect, lapply(abundance_list, colnames))
integrated_compositions <- do.call(rbind, lapply(abundance_list, `[`, common_columns))
rownames(integrated_compositions) <- gsub("\\.", "_", rownames(integrated_compositions))

integrated_compositions <- integrated_compositions[colnames(st_obj), common_columns, drop = FALSE]

write.csv(
  integrated_compositions,
  file = file.path(outdir, "cell2location_integrated_compositions.csv"),
  quote = FALSE
)

# ------------------------------------------------------------------------------
# 4. ILR transformation
# ------------------------------------------------------------------------------

base_ilr <- ilrBase(x = integrated_compositions, method = "basic")
cell_ilr <- as.matrix(ilr(integrated_compositions, base_ilr))
colnames(cell_ilr) <- paste0("ILR_", seq_len(ncol(cell_ilr)))

# ------------------------------------------------------------------------------
# 5. Multi-k Louvain clustering
# ------------------------------------------------------------------------------

cluster_info <- purrr::map(k_vect, function(k) {
  message("Building SNN graph for k = ", k)

  snn_graph <- scran::buildSNNGraph(
    x = t(as.matrix(as.data.frame(cell_ilr))),
    k = k
  )

  message("Running Louvain clustering for k = ", k)

  clust_louvain <- igraph::cluster_louvain(snn_graph, resolution = selected_resolution)

  tibble(
    spot_id = rownames(cell_ilr),
    cluster = clust_louvain$membership,
    k_label = paste0("k_", k)
  )
})

cluster_info <- bind_rows(cluster_info) %>%
  tidyr::pivot_wider(names_from = k_label, values_from = cluster)

write.csv(
  cluster_info,
  file = file.path(outdir, "niche_cluster_assignments.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. UMAP visualization
# ------------------------------------------------------------------------------

comp_umap <- uwot::umap(
  cell_ilr,
  n_neighbors = 30,
  n_epochs = 1000,
  metric = "cosine",
  min_dist = 0.1
) %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(row_id = rownames(cell_ilr))

write.csv(
  comp_umap,
  file = file.path(outdir, "niche_umap_coordinates.csv"),
  row.names = FALSE
)

plot_df <- comp_umap %>%
  left_join(cluster_info, by = c("row_id" = "spot_id"))

n_niches <- length(unique(plot_df[[paste0("k_", selected_k)]]))
plot_cols <- rep(defined_cols, length.out = n_niches)

p_umap <- ggplot(
  plot_df,
  aes(x = UMAP1, y = UMAP2, color = factor(.data[[paste0("k_", selected_k)]]))
) +
  ggrastr::geom_point_rast(size = 0.3, alpha = 0.8) +
  scale_color_manual(values = plot_cols) +
  theme_classic(base_size = 14) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Niche cluster",
    title = "Spatial niche UMAP"
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = file.path(outdir, "niche_umap.pdf"),
  plot = p_umap,
  width = 8,
  height = 6
)

# ------------------------------------------------------------------------------
# 7. Build niche metadata
# ------------------------------------------------------------------------------

selected_cluster_col <- paste0("k_", selected_k)

cluster_info_selected <- plot_df %>%
  dplyr::select(row_id, all_of(selected_cluster_col)) %>%
  dplyr::rename(niche = all_of(selected_cluster_col)) %>%
  mutate(ct_niche = paste0("niche_", niche))

# ------------------------------------------------------------------------------
# 8. Summarize niche composition
# ------------------------------------------------------------------------------

niche_summary_pat <- integrated_compositions %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id, values_to = "ct_prop", names_to = "cell_type") %>%
  left_join(cluster_info_selected, by = "row_id") %>%
  mutate(orig.ident = stringr::str_extract(row_id, "PT\\d+")) %>%
  group_by(orig.ident, ct_niche, cell_type) %>%
  summarize(median_ct_prop = median(ct_prop), .groups = "drop")

write.table(
  niche_summary_pat,
  file = file.path(outdir, "niche_summary_per_patient.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

niche_summary <- niche_summary_pat %>%
  group_by(ct_niche, cell_type) %>%
  summarize(patient_median_ct_prop = median(median_ct_prop), .groups = "drop")

# ------------------------------------------------------------------------------
# 9. Wilcoxon test for characteristic cell types
# ------------------------------------------------------------------------------

run_wilcox_up <- function(prop_data) {
  prop_data_group <- unique(prop_data$ct_niche)
  names(prop_data_group) <- prop_data_group

  purrr::map(prop_data_group, function(g) {
    test_data <- prop_data %>%
      mutate(
        test_group = ifelse(ct_niche == g, "target", "rest"),
        test_group = factor(test_group, levels = c("target", "rest"))
      )

    wilcox.test(median_ct_prop ~ test_group, data = test_data, alternative = "greater") %>%
      broom::tidy() %>%
      mutate(ct_niche = g)
  }) %>%
    bind_rows()
}

wilcoxon_res <- niche_summary_pat %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = purrr::map(data, run_wilcox_up)) %>%
  dplyr::select(cell_type, wres) %>%
  unnest(wres) %>%
  ungroup() %>%
  mutate(
    p_corr = p.adjust(p.value, method = "BH"),
    significant = ifelse(p_corr <= 0.05, "*", "")
  )

write.table(
  wilcoxon_res,
  file = file.path(outdir, "niche_wilcoxon_results.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# ------------------------------------------------------------------------------
# 10. Heatmap of niche-defining cell compositions
# ------------------------------------------------------------------------------

niche_order <- unique(niche_summary$ct_niche)
celltype_order_valid <- intersect(celltype_order, unique(niche_summary$cell_type))

p_heatmap <- niche_summary %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order_valid),
    ct_niche = factor(ct_niche, levels = rev(niche_order))
  ) %>%
  group_by(cell_type) %>%
  mutate(
    scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop)) / sd(patient_median_ct_prop)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = cell_type, y = ct_niche, fill = scaled_pat_median)) +
  geom_tile() +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = NULL, fill = "Scaled\nabundance")

ggsave(
  filename = file.path(outdir, "niche_composition_heatmap.pdf"),
  plot = p_heatmap,
  width = 10,
  height = 6
)

# ------------------------------------------------------------------------------
# 11. Niche proportions
# ------------------------------------------------------------------------------

cluster_counts <- cluster_info_selected %>%
  group_by(ct_niche) %>%
  summarise(nspots = dplyr::n(), .groups = "drop") %>%
  mutate(prop_spots = nspots / sum(nspots))

write.csv(
  cluster_counts,
  file = file.path(outdir, "niche_proportion_summary.csv"),
  row.names = FALSE
)

p_bar <- ggplot(cluster_counts, aes(x = ct_niche, y = prop_spots)) +
  geom_col() +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Proportion of spots", title = "Niche proportion summary")

ggsave(
  filename = file.path(outdir, "niche_proportion_barplot.pdf"),
  plot = p_bar,
  width = 6,
  height = 4
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_02_niche_analysis.txt")
)