#!/usr/bin/env Rscript
# 00b_dimred_clustering.R
# Dimensionality reduction (PCA, UMAP, t-SNE) and clustering
# HIGH-QUALITY PLOTTING VERSION

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
})

source("scripts/spata2/plot_themes.R")

OBJECTS_DIR <- "analysis/objects"
PLOTS_DIR <- "analysis/plots/dimred"
LOGS_DIR <- "analysis/logs"

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGS_DIR, paste0("00b_dimred_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Dimensionality Reduction and Clustering ===\n")

SAMPLE_NAMES <- c("D1_27038840_HCCboneMet_seg", "SRout_27038840_HCC_tumour")

for (sample_name in SAMPLE_NAMES) {
  cat("\n--- Processing:", sample_name, "---\n")
  
  norm_file <- file.path(OBJECTS_DIR, paste0(sample_name, "_norm_lognorm.rds"))
  
  if (!file.exists(norm_file)) {
    cat("  Normalized object not found, skipping.\n")
    next
  }
  
  spata_obj <- readRDS(norm_file)
  cat("  Loaded:", nObs(spata_obj), "spots\n")
  
  # Find HVGs
  cat("  Finding highly variable genes...\n")
  spata_obj <- tryCatch({
    findHVGs(spata_obj, n_hvgs = 2000)
  }, error = function(e) {
    cat("    HVG detection failed:", e$message, "\n")
    spata_obj
  })
  
  # PCA
  cat("  Computing PCA (30 components)...\n")
  spata_obj <- tryCatch({
    runPCA(spata_obj, n_pcs = 30)
  }, error = function(e) {
    cat("    PCA failed:", e$message, "\n")
    spata_obj
  })
  
  # UMAP
  cat("  Computing UMAP...\n")
  spata_obj <- tryCatch({
    runUMAP(spata_obj, n_pcs = 20)
  }, error = function(e) {
    cat("    UMAP failed:", e$message, "\n")
    spata_obj
  })
  
  # t-SNE
  cat("  Computing t-SNE...\n")
  spata_obj <- tryCatch({
    runTSNE(spata_obj, n_pcs = 20)
  }, error = function(e) {
    cat("    t-SNE failed:", e$message, "\n")
    spata_obj
  })
  
  # K-means clustering
  cat("  K-means clustering (k=7)...\n")
  spata_obj <- tryCatch({
    # Get PCA coordinates
    pca_df <- getDimRedDf(spata_obj, method = "pca")
    
    if (!is.null(pca_df) && nrow(pca_df) > 0) {
      pca_cols <- grep("^PC", colnames(pca_df), value = TRUE)[1:min(20, sum(grepl("^PC", colnames(pca_df))))]
      pca_mat <- as.matrix(pca_df[, pca_cols])
      
      set.seed(42)
      kmeans_result <- kmeans(pca_mat, centers = 7, nstart = 25, iter.max = 100)
      
      cluster_df <- data.frame(
        barcodes = pca_df$barcodes,
        kmeans_cluster = paste0("Cluster_", kmeans_result$cluster)
      )
      
      spata_obj <- addMetaData(spata_obj, meta_df = cluster_df)
      cat("    ✓ K-means clustering complete\n")
    }
    spata_obj
  }, error = function(e) {
    cat("    K-means failed:", e$message, "\n")
    spata_obj
  })
  
  # Generate plots
  cat("  Generating dimension reduction plots...\n")
  
  sample_plot_dir <- file.path(PLOTS_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    # Get UMAP coordinates
    umap_df <- getDimRedDf(spata_obj, method = "umap")
    tsne_df <- getDimRedDf(spata_obj, method = "tsne")
    
    meta <- getMetaDf(spata_obj)
    
    if (!is.null(umap_df) && "kmeans_cluster" %in% colnames(meta)) {
      umap_df <- merge(umap_df, meta[, c("barcodes", "kmeans_cluster")], by = "barcodes")
      
      p_umap <- ggplot(umap_df, aes(x = umap1, y = umap2, color = kmeans_cluster)) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_manual(values = DISCRETE_COLORS[1:7]) +
        HQ_THEME() +
        labs(title = paste(sample_name, "- UMAP"), color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      
      save_plot(p_umap, file.path(sample_plot_dir, "umap_clusters.png"), width = 10, height = 8)
    }
    
    if (!is.null(tsne_df) && "kmeans_cluster" %in% colnames(meta)) {
      tsne_df <- merge(tsne_df, meta[, c("barcodes", "kmeans_cluster")], by = "barcodes")
      
      p_tsne <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = kmeans_cluster)) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_manual(values = DISCRETE_COLORS[1:7]) +
        HQ_THEME() +
        labs(title = paste(sample_name, "- t-SNE"), color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      
      save_plot(p_tsne, file.path(sample_plot_dir, "tsne_clusters.png"), width = 10, height = 8)
    }
    
    cat("    ✓ Plots saved\n")
    
  }, error = function(e) {
    cat("    Plot generation failed:", e$message, "\n")
  })
  
  # Save updated object
  saveRDS(spata_obj, norm_file)
  cat("  Saved updated object\n")
  
  rm(spata_obj)
  gc()
}

cat("\n=== Dimensionality Reduction Complete ===\n")
sink()
