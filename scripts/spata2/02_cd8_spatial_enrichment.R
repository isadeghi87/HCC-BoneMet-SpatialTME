#!/usr/bin/env Rscript
# 02_cd8_spatial_enrichment.R
# CD8 spatial enrichment analysis - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

PLOT_DIR <- "analysis/plots/cd8_spatial"
TABLE_DIR <- "analysis/tables/cd8_spatial"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("02_cd8_spatial_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== CD8 Spatial Enrichment Analysis ===\n")

CD8_MARKERS <- c("CD8A", "CD8B")

analyze_cd8_spatial <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  output_prefix <- get_output_prefix(sample_name)
  tissue_label <- get_tissue_label(get_tissue_from_visium(sample_name))
  cat("Output prefix:", output_prefix, "\n")
  obj <- readRDS(obj_path)
  cat("Loaded object with", nrow(obj@meta_obs), "spots\n")
  
  available_genes <- rownames(getMatrix(obj))
  cd8_available <- CD8_MARKERS[CD8_MARKERS %in% available_genes]
  cd8_missing <- CD8_MARKERS[!CD8_MARKERS %in% available_genes]
  
  cat("CD8 markers available:", paste(cd8_available, collapse = ", "), "\n")
  if (length(cd8_missing) > 0) cat("CD8 markers missing:", paste(cd8_missing, collapse = ", "), "\n")
  
  if (length(cd8_available) == 0) {
    cat("ERROR: No CD8 markers found. Skipping.\n")
    return(NULL)
  }
  
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  cd8_expr <- expr_mtr[cd8_available, , drop = FALSE]
  cd8_score <- colMeans(as.matrix(cd8_expr))
  
  obj@meta_obs$cd8_score <- cd8_score[obj@meta_obs$barcodes]
  for (marker in cd8_available) {
    obj@meta_obs[[marker]] <- as.numeric(expr_mtr[marker, obj@meta_obs$barcodes])
  }
  
  q75 <- quantile(cd8_score, 0.75, na.rm = TRUE)
  obj@meta_obs$cd8_high <- ifelse(obj@meta_obs$cd8_score >= q75, "CD8-high", "CD8-low")
  
  cat("CD8 score quantiles:\n")
  print(quantile(cd8_score, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  cat("CD8-high spots:", sum(obj@meta_obs$cd8_high == "CD8-high"), "\n")
  cat("CD8-low spots:", sum(obj@meta_obs$cd8_high == "CD8-low"), "\n")
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # CD8 score spatial map
  p_score <- ggplot(plot_df, aes(x = x, y = y, color = cd8_score)) +
    geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
    scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
    coord_fixed() + theme_publication() +
    theme(legend.key.width = unit(1.5, "cm")) +
    labs(title = paste(tissue_label, "- CD8 Score"), color = "CD8 Score")
  save_plot(p_score, file.path(sample_plot_dir, paste0(output_prefix, "_ST_cd8_score_spatial.png")), width = 12, height = 10)
  
  # CD8-high/low map
  p_high <- ggplot(plot_df, aes(x = x, y = y, color = cd8_high)) +
    geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
    scale_color_manual(values = c("CD8-high" = "#d62728", "CD8-low" = "#d9d9d9")) +
    coord_fixed() + theme_publication() +
    labs(title = paste(tissue_label, "- CD8-high Regions"), color = "Status") +
    guides(color = guide_legend(override.aes = list(size = 4)))
  save_plot(p_high, file.path(sample_plot_dir, paste0(output_prefix, "_ST_cd8_high_regions.png")), width = 12, height = 10)
  
  # Individual markers
  for (marker in cd8_available) {
    p_marker <- ggplot(plot_df, aes(x = x, y = y, color = .data[[marker]])) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
      coord_fixed() + theme_publication() +
      labs(title = paste(tissue_label, "-", marker), color = marker)
    save_plot(p_marker, file.path(sample_plot_dir, paste0(output_prefix, "_ST_", marker, "_spatial.png")), width = 12, height = 10)
  }
  
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    cd8_mean = mean(obj@meta_obs$cd8_score, na.rm = TRUE),
    cd8_median = median(obj@meta_obs$cd8_score, na.rm = TRUE),
    cd8_q75 = q75,
    cd8_high_count = sum(obj@meta_obs$cd8_high == "CD8-high"),
    cd8_low_count = sum(obj@meta_obs$cd8_high == "CD8-low")
  )
  write.csv(stats_df, file.path(sample_table_dir, paste0(output_prefix, "_ST_cd8_spatial_stats.csv")), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Updated object saved.\n")
  return(stats_df)
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

all_stats <- list()
for (s in samples) {
  if (file.exists(s$path)) {
    stats <- tryCatch(analyze_cd8_spatial(s$path, s$name), error = function(e) { cat("Error:", e$message, "\n"); NULL })
    if (!is.null(stats)) all_stats[[s$name]] <- stats
  }
}

if (length(all_stats) > 0) {
  combined_stats <- do.call(rbind, all_stats)
  write.csv(combined_stats, file.path(TABLE_DIR, "cd8_spatial_stats_combined.csv"), row.names = FALSE)
}

cat("\n=== CD8 Spatial Enrichment Complete ===\n")
sink()
