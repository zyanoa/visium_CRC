# ==============================================================================
# Script: 05_misty_group_visualization.R
# Project: Spatial transcriptomic atlas of dMMR and pMMR colorectal cancer
# Purpose:
#   1. Collect MISTy results across all sections
#   2. Summarize view contributions and explained variance (R2)
#   3. Compare dMMR and pMMR groups
#   4. Generate group-level interaction heatmaps
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mistyR)
  library(broom)
  library(readr)
  library(tidyr)
  library(dplyr)
})

source("00_utils/misty_utilities.R")

set.seed(1234)

CONFIG_FILE <- Sys.getenv("CRC_MMR_CONFIG", unset = "00_utils/project_config.R")
source(CONFIG_FILE)
load_project_config_message()

# ------------------------------------------------------------------------------
# 1. Input and output
# ------------------------------------------------------------------------------

misty_out_folder <- file.path(RESULTS_DIR, "03_spatial_ecology", "misty")
outdir <- file.path(RESULTS_DIR, "03_spatial_ecology", "misty_group_summary")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Sample metadata
# ------------------------------------------------------------------------------

dmmr_patients <- c("PT18", "PT32", "PT9")
pmmr_patients <- c("PT34", "PT55", "PT57")

cell_order <- list(
  predictor = c(
    "TCD4", "NK", "TCD8", "Mast", "TZBTB16", "Endo", "Schwann",
    "Fibro", "Peri", "SmoothMuscle", "Mono", "Granulo", "Macro",
    "B", "Epi", "ILC", "γδ_T", "DC", "Plasma"
  ),
  target = c(
    "TZBTB16", "B", "TCD4", "TCD8", "Mast", "NK", "DC",
    "γδ_T", "Plasma", "Epi", "ILC", "Mono", "Granulo", "Macro",
    "Schwann", "Endo", "SmoothMuscle", "Fibro", "Peri"
  )
)

# ------------------------------------------------------------------------------
# 3. Collect MISTy output directories
# ------------------------------------------------------------------------------

sample_dirs <- list.dirs(misty_out_folder, full.names = TRUE, recursive = FALSE)

misty_result_dirs <- unlist(lapply(sample_dirs, function(sample_dir) {
  subdirs <- list.dirs(sample_dir, full.names = TRUE, recursive = FALSE)
  subdirs[grepl("_cell2location$", basename(subdirs))]
}))

if (length(misty_result_dirs) == 0) {
  stop("No MISTy result directories found.")
}

misty_res <- collect_results(misty_result_dirs)

# ------------------------------------------------------------------------------
# 4. Prepare metadata for improvements and importances
# ------------------------------------------------------------------------------

sample_importances <- misty_res$importances
sample_importances$sample <- gsub(".*/(PT[0-9]+_[0-9]+)/.*", "\\1", sample_importances$sample)
sample_importances$patient <- gsub("^(PT[0-9]+)_.*", "\\1", sample_importances$sample)
sample_importances$group <- case_when(
  sample_importances$patient %in% dmmr_patients ~ "dMMR",
  sample_importances$patient %in% pmmr_patients ~ "pMMR",
  TRUE ~ NA_character_
)

R2_data <- misty_res$improvements %>%
  dplyr::filter(measure == "multi.R2")

R2_data$sample <- gsub(".*/(PT[0-9]+_[0-9]+)/.*", "\\1", R2_data$sample)
R2_data$patient <- gsub("^(PT[0-9]+)_.*", "\\1", R2_data$sample)
R2_data$group <- case_when(
  R2_data$patient %in% dmmr_patients ~ "dMMR",
  R2_data$patient %in% pmmr_patients ~ "pMMR",
  TRUE ~ NA_character_
)
colnames(R2_data)[colnames(R2_data) == "value"] <- "R2"

write.csv(R2_data, file = file.path(outdir, "misty_R2_all_samples.csv"), row.names = FALSE)
write.csv(sample_importances, file = file.path(outdir, "misty_importances_all_samples.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Quick global visualization
# ------------------------------------------------------------------------------

pdf(file.path(outdir, "view_contributions_all_samples.pdf"), width = 8, height = 6)
plot_view_contributions(misty_res)
dev.off()

pdf(file.path(outdir, "interaction_communities_all_samples_intra.pdf"), width = 8, height = 6)
plot_interaction_communities(misty_res, "intra")
dev.off()

# ------------------------------------------------------------------------------
# 6. Helper functions
# ------------------------------------------------------------------------------

plot_r2_box <- function(df, target_levels = NULL) {
  if (!is.null(target_levels)) {
    df$target <- factor(df$target, levels = target_levels)
  }

  ggplot(df, aes(x = target, y = R2)) +
    geom_boxplot() +
    geom_point(aes(color = group), position = position_jitter(width = 0.15)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    ) +
    ylab("Explained variance (R2)") +
    xlab(NULL)
}

plot_r2_tile <- function(df, target_levels = NULL) {
  if (!is.null(target_levels)) {
    df$target <- factor(df$target, levels = target_levels)
  }

  ggplot(df, aes(x = target, y = sample, fill = R2)) +
    geom_tile() +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    ) +
    scale_fill_gradient(low = "black", high = "yellow") +
    xlab(NULL) +
    ylab(NULL)
}

plot_importance_heatmap <- function(view_importance, cell_order, color_fill = "#8DA0CB", title = NULL) {
  p <- view_importance %>%
    mutate(
      Predictor = factor(Predictor, levels = cell_order$predictor),
      Target = factor(Target, levels = cell_order$target)
    ) %>%
    ggplot(aes(x = Target, y = Predictor, fill = median_importance)) +
    geom_tile() +
    scale_fill_gradient2(high = color_fill, midpoint = 0, low = "white", na.value = "grey") +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    ) +
    xlab(NULL) +
    ylab(NULL)

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  p
}

# ------------------------------------------------------------------------------
# 7. Group-wise summarization
# ------------------------------------------------------------------------------

groups <- c("dMMR", "pMMR")
R2_threshold <- 50

for (grp in groups) {
  message("Processing group: ", grp)

  current_R2_data <- subset(R2_data, group == grp)
  current_importances <- subset(sample_importances, group == grp)

  group_outdir <- file.path(outdir, grp)
  if (!dir.exists(group_outdir)) dir.create(group_outdir, recursive = TRUE)

  write.csv(current_R2_data, file = file.path(group_outdir, "R2_results.csv"), row.names = FALSE)
  write.csv(current_importances, file = file.path(group_outdir, "sample_importances_raw.csv"), row.names = FALSE)

  p_r2_box <- plot_r2_box(current_R2_data, target_levels = cell_order$target) + ggtitle(grp)
  p_r2_tile <- plot_r2_tile(current_R2_data, target_levels = cell_order$target) + ggtitle(grp)

  ggsave(file.path(group_outdir, "R2_boxplot.pdf"), p_r2_box, width = 6, height = 4)
  ggsave(file.path(group_outdir, "R2_tileplot.pdf"), p_r2_tile, width = 7, height = 7)

  best_performers <- current_R2_data %>%
    dplyr::select(target, sample, R2) %>%
    dplyr::filter(R2 >= R2_threshold) %>%
    mutate(best_performer = TRUE)

  importances_filt <- current_importances %>%
    left_join(
      best_performers %>% dplyr::select(-R2),
      by = c("Target" = "target", "sample")
    ) %>%
    dplyr::filter(!is.na(best_performer))

  write.csv(importances_filt, file = file.path(group_outdir, "sample_importances_filtered.csv"), row.names = FALSE)

  summarized_interactions <- importances_filt %>%
    group_by(view, Predictor, Target) %>%
    summarize(
      median_importance = median(Importance, na.rm = TRUE),
      .groups = "drop"
    )

  importance_test <- importances_filt %>%
    na.omit() %>%
    dplyr::select(view, Predictor, Target, Importance) %>%
    group_by(view, Predictor, Target) %>%
    nest() %>%
    mutate(
      wres = purrr::map(data, function(dat) {
        wilcox.test(dat$Importance, mu = 0, alternative = "greater") %>%
          broom::tidy()
      })
    ) %>%
    dplyr::select(-data) %>%
    unnest(wres)

  summarized_interactions <- summarized_interactions %>%
    left_join(importance_test, by = c("view", "Predictor", "Target"))

  write.csv(
    summarized_interactions,
    file = file.path(group_outdir, "summarized_interactions.csv"),
    row.names = FALSE
  )

  intra_df <- summarized_interactions %>% filter(view == "intra")
  juxta_df <- summarized_interactions %>% filter(view == "juxta_2")
  para_df <- summarized_interactions %>% filter(view == "para_5")

  p_intra <- plot_importance_heatmap(intra_df, cell_order, color_fill = "darkgreen", title = paste0(grp, "_Intra"))
  p_juxta <- plot_importance_heatmap(juxta_df, cell_order, color_fill = "orange", title = paste0(grp, "_Juxta_2"))
  p_para <- plot_importance_heatmap(para_df, cell_order, color_fill = "#8DA0CB", title = paste0(grp, "_Para_5"))

  ggsave(file.path(group_outdir, "misty_importances_intra.pdf"), p_intra, width = 7, height = 7)
  ggsave(file.path(group_outdir, "misty_importances_juxta.pdf"), p_juxta, width = 7, height = 7)
  ggsave(file.path(group_outdir, "misty_importances_para.pdf"), p_para, width = 7, height = 7)
}

# ------------------------------------------------------------------------------
# 8. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "sessionInfo_05_misty_group_visualization.txt")
)