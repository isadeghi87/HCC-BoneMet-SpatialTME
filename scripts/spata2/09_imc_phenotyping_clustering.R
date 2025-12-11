#!/usr/bin/env Rscript
# 09_imc_phenotyping_clustering.R
# IMC FlowSOM clustering and cell phenotyping - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(FlowSOM)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

OBJ_DIR <- "analysis/objects/imc"
PLOT_DIR <- "analysis/plots/imc"
TABLE_DIR <- "analysis/tables/imc"
LOG_DIR <- "analysis/logs"

log_file <- file.path(LOG_DIR, paste0("09_imc_clustering_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== IMC FlowSOM Clustering ===\n")

IMMUNE_MARKERS <- c("CD45", "CD3", "CD4", "CD8", "CD20", "CD68", "CD163", "CD56", "CD14")
FUNCTIONAL_MARKERS <- c("Ki67", "PD1", "PDL1", "CTLA4", "LAG3", "TIM3", "GZMB", "PRF1")

process_imc_clustering <- function(sample_path) {
  sample_name <- gsub("_processed\\.rds$", "", basename(sample_path))
  cat("\n--- Processing:", sample_name, "---\n")
  
  # Extract standardized naming components
  # Parse sample name to get base (e.g., "27038840_Biopsy_Spine" from "27038840_Biopsy_Spine_1_1")
  base_sample <- sub("_\\d+_\\d+$", "", sample_name)
  roi_match <- regmatches(sample_name, regexpr("\\d+_\\d+$", sample_name))
  roi_num <- if (length(roi_match) > 0) as.integer(sub("_.*", "", roi_match)) else 1
  roi_label <- paste0("ROI", roi_num)
  
  # Get standardized output prefix and tissue label
  output_prefix <- get_output_prefix(base_sample)
  tissue_type <- get_tissue_from_imc(base_sample)
  tissue_label <- if (!is.na(tissue_type)) get_tissue_label(tissue_type) else sample_name
  
  cat("  Base sample:", base_sample, "\n")
  cat("  ROI:", roi_label, "\n")
  cat("  Output prefix:", output_prefix, "\n")
  cat("  Tissue label:", tissue_label, "\n")
  
  dat <- readRDS(sample_path)
  expr_norm <- dat$expr_norm
  coords <- dat$coords
  
  cat("Cells:", nrow(expr_norm), ", Markers:", ncol(expr_norm), "\n")
  
  # Subsample for large datasets
  max_cells <- 100000
  if (nrow(expr_norm) > max_cells) {
    cat("Subsampling to", max_cells, "cells for clustering\n")
    idx <- sample(nrow(expr_norm), max_cells)
    expr_sub <- expr_norm[idx, ]
    coords_sub <- if (!is.null(coords)) coords[idx, ] else NULL
  } else {
    expr_sub <- expr_norm
    coords_sub <- coords
  }
  
  markers <- colnames(expr_sub)
  clustering_markers <- markers[!markers %in% c("Ki67", "PD1", "PDL1", "CTLA4", "LAG3", "TIM3")]
  cat("Clustering markers:", length(clustering_markers), "\n")
  
  cat("Running FlowSOM...\n")
  fsom <- tryCatch({
    FlowSOM::BuildSOM(
      data = expr_sub[, clustering_markers, drop = FALSE],
      colsToUse = seq_along(clustering_markers),
      xdim = 10, ydim = 10, rlen = 10
    )
  }, error = function(e) { cat("FlowSOM error:", e$message, "\n"); NULL })
  
  if (is.null(fsom)) {
    cat("Falling back to K-means...\n")
    set.seed(42)
    km <- kmeans(expr_sub[, clustering_markers], centers = 15, nstart = 10, iter.max = 100)
    clusters <- km$cluster
  } else {
    cat("Performing metaclustering...\n")
    mc <- FlowSOM::MetaClustering(fsom$map$codes, method = "metaClustering_consensus", max = 15)
    clusters <- mc[fsom$map$mapping[, 1]]
  }
  
  dat$clusters_sub <- clusters
  dat$subsample_idx <- if (nrow(expr_norm) > max_cells) idx else seq_len(nrow(expr_norm))
  
  cat("Clusters:", length(unique(clusters)), "\n")
  print(table(clusters))
  
  # Cluster profiles
  cluster_profiles <- aggregate(expr_sub, by = list(cluster = clusters), FUN = mean)
  rownames(cluster_profiles) <- paste0("C", cluster_profiles$cluster)
  cluster_profiles$cluster <- NULL
  
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # HIGH-QUALITY heatmap
  heatmap_file <- paste0(output_prefix, "_IMC_", roi_label, "_cluster_heatmap.png")
  png(file.path(sample_plot_dir, heatmap_file), width = 1400, height = 1000, res = 150)
  tryCatch({
    pheatmap(
      as.matrix(cluster_profiles),
      scale = "column",
      cluster_rows = TRUE, cluster_cols = TRUE,
      color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
      main = paste(tissue_label, "-", roi_label, "- Cluster Profiles"),
      fontsize = 10, fontsize_row = 9, fontsize_col = 8
    )
  }, error = function(e) cat("Heatmap error:", e$message, "\n"))
  dev.off()
  
  # Spatial cluster plot
  if (!is.null(coords_sub)) {
    plot_df <- data.frame(coords_sub, cluster = factor(clusters))
    
    p_spatial <- ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.15, alpha = 0.7, shape = 16, stroke = 0) +
      scale_color_manual(values = c(DISCRETE_COLORS, brewer.pal(8, "Dark2"))) +
      coord_fixed() + theme_publication() +
      labs(title = paste(tissue_label, "-", roi_label, "- Cluster Spatial"), color = "Cluster") +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
    spatial_file <- paste0(output_prefix, "_IMC_", roi_label, "_cluster_spatial.png")
    save_plot(p_spatial, file.path(sample_plot_dir, spatial_file), width = 12, height = 10)
  }
  
  # Basic annotation
  cat("\nAnnotating clusters...\n")
  annotations <- sapply(seq_len(nrow(cluster_profiles)), function(i) {
    profile <- as.numeric(cluster_profiles[i, ])
    names(profile) <- colnames(cluster_profiles)
    
    cd45_cols <- grep("CD45", names(profile), value = TRUE, ignore.case = TRUE)
    cd3_cols <- grep("CD3", names(profile), value = TRUE, ignore.case = TRUE)
    cd8_cols <- grep("CD8", names(profile), value = TRUE, ignore.case = TRUE)
    cd68_cols <- grep("CD68", names(profile), value = TRUE, ignore.case = TRUE)
    
    cd45_expr <- if (length(cd45_cols) > 0) mean(profile[cd45_cols]) else 0
    cd3_expr <- if (length(cd3_cols) > 0) mean(profile[cd3_cols]) else 0
    cd8_expr <- if (length(cd8_cols) > 0) mean(profile[cd8_cols]) else 0
    cd68_expr <- if (length(cd68_cols) > 0) mean(profile[cd68_cols]) else 0
    
    if (cd45_expr > 1 && cd3_expr > 1 && cd8_expr > 1) return("CD8 T cell")
    if (cd45_expr > 1 && cd3_expr > 1) return("T cell")
    if (cd45_expr > 1 && cd68_expr > 1) return("Macrophage")
    if (cd45_expr > 1) return("Immune")
    return("Other")
  })
  
  dat$annotations <- data.frame(cluster = rownames(cluster_profiles), annotation = annotations)
  cat("Annotations:\n")
  print(table(annotations))
  
  profiles_file <- paste0(output_prefix, "_IMC_", roi_label, "_cluster_profiles.csv")
  annotations_file <- paste0(output_prefix, "_IMC_", roi_label, "_cluster_annotations.csv")
  write.csv(cluster_profiles, file.path(sample_table_dir, profiles_file))
  write.csv(dat$annotations, file.path(sample_table_dir, annotations_file), row.names = FALSE)
  
  saveRDS(dat, sample_path)
  cat("Done:", sample_name, "(saved as", output_prefix, "_IMC_", roi_label, ")\n")
}

processed_files <- list.files(OBJ_DIR, pattern = "_processed\\.rds$", full.names = TRUE)
cat("Found", length(processed_files), "processed files\n")

for (f in processed_files) {
  tryCatch(process_imc_clustering(f), error = function(e) cat("Error processing", f, ":", e$message, "\n"))
}

cat("\n=== Clustering Complete ===\n")
sink()
