# ==============================================================================
# File: MIA_function.R
# Purpose:
#   Utility function for multimodal intersection analysis (MIA) between
#   spatial transcriptomic regions and scRNA-seq cell types.
# ==============================================================================

zhao_MIA <- function(
    sp.diff,
    sc.diff,
    sample,
    outdir,
    sorted_regions = c(
      "dMMR_Cancer1_CEACAM5",
      "pMMR_Cancer1_CEACAM5",
      "dMMR_Cancer2_KRT8",
      "pMMR_Cancer2_KRT8",
      "dMMR_Normal_PIGR",
      "pMMR_Normal_PIGR",
      "dMMR_Stroma1_IGLC1",
      "pMMR_Stroma1_IGLC1",
      "dMMR_Stroma2_MYL9+ACTG2",
      "pMMR_Stroma2_MYL9+ACTG2",
      "dMMR_Stroma3_APOE",
      "pMMR_Stroma3_APOE",
      "dMMR_Stroma4_IGHM",
      "pMMR_Stroma4_IGHM",
      "dMMR_Stroma5_CXCL8",
      "pMMR_Stroma5_CXCL8"
    ),
    celltype_order = c(
      "Epi", "TCD4", "TCD8", "Tgd", "NK", "Mast", "TZBTB16",
      "B", "Plasma", "ILC", "Peri", "Endo", "Fibro",
      "SmoothMuscle", "Granulo", "Mono", "Macro", "DC", "Schwann"
    ),
    pvalue_cutoff = 1,
    qvalue_cutoff = 1,
    min_gssize = 1,
    max_gssize = 10000,
    padjust_method = "BH"
) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(foreach)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(grid)
  })

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  stopifnot(all(c("region", "gene") %in% colnames(sp.diff)))
  stopifnot(all(c("celltype", "gene") %in% colnames(sc.diff)))

  term2name <- data.frame(
    term = sp.diff$region,
    name = sp.diff$region,
    stringsAsFactors = FALSE
  )

  term2gene <- data.frame(
    sp.diff,
    stringsAsFactors = FALSE
  )

  cluster_list <- unique(sc.diff$celltype)
  enrich_list <- list()

  for (each_cluster in cluster_list) {
    diff_gene <- sc.diff %>%
      filter(celltype == each_cluster) %>%
      pull(gene) %>%
      as.character() %>%
      unique()

    enrich <- tryCatch(
      enricher(
        gene = diff_gene,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = padjust_method,
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
        minGSSize = min_gssize,
        maxGSSize = max_gssize
      ),
      error = function(e) NULL
    )

    if (!is.null(enrich) && nrow(enrich@result) > 0) {
      enrich_res <- enrich@result[, c("ID", "Description", "p.adjust")]
      enrich_res$sc.celltype <- each_cluster
      enrich_res$ES <- -log10(enrich_res$p.adjust)
      rownames(enrich_res) <- enrich_res$ID
      enrich_list[[each_cluster]] <- enrich_res
    }
  }

  if (length(enrich_list) == 0) {
    warning("No valid enrichment results were generated.")
    return(NULL)
  }

  combined_result <- do.call(rbind, enrich_list)
  write.csv(
    combined_result,
    file = file.path(outdir, paste0(sample, "_MIA.Result.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  all_terms <- NULL
  for (res in enrich_list) {
    idx <- which(!rownames(res) %in% names(all_terms))
    new_terms <- res$Description[idx]
    names(new_terms) <- rownames(res)[idx]
    all_terms <- c(all_terms, new_terms)
  }

  padj_df <- foreach(res = enrich_list, .combine = rbind) %do% {
    padj <- res[names(all_terms), "p.adjust"]
    padj[is.na(padj)] <- 1
    names(padj) <- names(all_terms)
    return(padj)
  }

  rownames(padj_df) <- names(enrich_list)
  padj_df <- t(as.matrix(padj_df))

  min_qadj_terms <- apply(padj_df, 1, min)
  min_qadj_terms <- sort(min_qadj_terms)

  selected_idx <- which(min_qadj_terms < 0.05 & seq_along(min_qadj_terms) <= 20)
  if (length(selected_idx) == 0) {
    selected_idx <- seq_len(min(20, length(min_qadj_terms)))
  }

  plot_terms <- names(min_qadj_terms)[selected_idx]

  if (length(plot_terms) <= 1) {
    warning("Not enough enriched terms to generate heatmap.")
    return(
      list(
        result_table = combined_result,
        padj_matrix = padj_df,
        plot_matrix = NULL
      )
    )
  }

  plot_data <- padj_df[plot_terms, , drop = FALSE]
  plot_data <- -log10(plot_data)
  plot_data[plot_data > 8] <- 8
  plot_data <- as.data.frame(plot_data)

  valid_cols <- intersect(celltype_order, colnames(plot_data))
  valid_rows <- intersect(sorted_regions, rownames(plot_data))

  plot_data <- plot_data[valid_rows, valid_cols, drop = FALSE]

  split_by_region <- gsub("^(dMMR_|pMMR_)", "", rownames(plot_data))
  split_by_region <- factor(split_by_region, levels = unique(split_by_region))

  colors_vals <- colorRampPalette(c("#FFF4ED", "darkred"))

  cell_size <- unit(1, "cm")
  heatmap_width <- ncol(plot_data) * cell_size
  heatmap_height <- nrow(plot_data) * cell_size

  ht <- Heatmap(
    as.matrix(plot_data),
    col = colorRamp2(breaks = 0:8, colors = colors_vals(9)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "Cell type",
    column_title_side = "bottom",
    rect_gp = gpar(col = "black"),
    show_row_names = TRUE,
    row_split = split_by_region,
    width = heatmap_width,
    height = heatmap_height,
    row_names_gp = gpar(fontsize = 16),
    column_names_gp = gpar(fontsize = 16),
    heatmap_legend_param = list(
      title = "Adjusted P\n(-log10)",
      title_position = "leftcenter",
      at = c(0, 2, 4, 6, 8),
      labels = c("0", "2", "4", "6", ">8"),
      legend_width = unit(0.5, "npc"),
      legend_direction = "horizontal",
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 16)
    )
  )

  png(
    filename = file.path(outdir, paste0(sample, ".cluster.MIA.enrich.heatmap.png")),
    type = "cairo-png",
    width = 13 * 200,
    height = 15 * 200,
    res = 200
  )
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()

  pdf(
    file = file.path(outdir, paste0(sample, ".cluster.MIA.enrich.heatmap.pdf")),
    width = 13,
    height = 15
  )
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()

  return(
    list(
      result_table = combined_result,
      padj_matrix = padj_df,
      plot_matrix = plot_data
    )
  )
}