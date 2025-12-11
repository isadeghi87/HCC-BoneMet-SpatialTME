#!/usr/bin/env Rscript
# 05_salie_immunotype_markers.R
# Salié Immunotype classification and TLS analysis - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")

PLOT_DIR <- "analysis/plots/salie"
TABLE_DIR <- "analysis/tables/salie"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("05_salie_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Salié Immunotype Analysis ===\n")

COL4A1_MARKERS <- c("COL4A1", "COL4A2")
DC_MARKERS <- c("CD1C", "CLEC9A", "CD141", "XCR1", "LAMP3")
TLS_MARKERS <- c("CXCL13", "CCL19", "CCL21", "CD79A", "MS4A1", "CD20")
EXCLUSION_MARKERS <- c("CTNNB1", "TGFB1", "COL4A1")

calc_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

analyze_salie <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  cat("COL4A1 markers:", paste(COL4A1_MARKERS[COL4A1_MARKERS %in% genes], collapse = ","), "\n")
  cat("DC markers:", paste(DC_MARKERS[DC_MARKERS %in% genes], collapse = ","), "\n")
  cat("TLS markers:", paste(TLS_MARKERS[TLS_MARKERS %in% genes], collapse = ","), "\n")
  
  obj@meta_obs$col4a1_score <- calc_score(expr_mtr, COL4A1_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$dc_score <- calc_score(expr_mtr, DC_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$tls_score <- calc_score(expr_mtr, TLS_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$exclusion_score <- calc_score(expr_mtr, EXCLUSION_MARKERS, genes)[obj@meta_obs$barcodes]
  
  if (!all(is.na(obj@meta_obs$tls_score))) {
    tls_thresh <- quantile(obj@meta_obs$tls_score, 0.9, na.rm = TRUE)
    obj@meta_obs$tls_hotspot <- ifelse(obj@meta_obs$tls_score >= tls_thresh, "TLS-high", "TLS-low")
    cat("TLS hotspots (>90th percentile):", sum(obj@meta_obs$tls_hotspot == "TLS-high"), "\n")
  }
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  scores <- c("col4a1_score", "dc_score", "tls_score", "exclusion_score")
  titles <- c("COL4A1 (Exclusion)", "Dendritic Cells", "TLS Signature", "Exclusion Score")
  
  plots <- list()
  for (i in seq_along(scores)) {
    if (!all(is.na(plot_df[[scores[i]]]))) {
      plots[[i]] <- ggplot(plot_df, aes(x = x, y = y, color = .data[[scores[i]]])) +
        geom_point(size = 0.25, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm")) +
        labs(title = titles[i], color = "Score")
    }
  }
  
  if (length(plots) > 0) {
    combined <- wrap_plots(plots[!sapply(plots, is.null)], ncol = 2) +
      plot_annotation(title = paste(sample_name, "- Salié Immunotype"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, "salie_combined.png"), width = 14, height = 12)
  }
  
  # TLS hotspot map
  if ("tls_hotspot" %in% names(plot_df)) {
    p_tls <- ggplot(plot_df, aes(x = x, y = y, color = tls_hotspot)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_manual(values = c("TLS-high" = "#d62728", "TLS-low" = "#d9d9d9")) +
      coord_fixed() + theme_publication() +
      labs(title = paste(sample_name, "- TLS Hotspots (Top 10%)"), color = "Status") +
      guides(color = guide_legend(override.aes = list(size = 4)))
    save_plot(p_tls, file.path(sample_plot_dir, "tls_hotspots.png"), width = 12, height = 10)
  }
  
  # Co-localisation: TLS with CD8-high
  if ("cd8_high" %in% names(obj@meta_obs) && "tls_hotspot" %in% names(obj@meta_obs)) {
    cat("TLS-CD8 co-localisation:\n")
    contingency <- table(obj@meta_obs$tls_hotspot, obj@meta_obs$cd8_high)
    print(contingency)
  }
  
  sample_table_dir <- file.path(TABLE_DIR, sample_name)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    col4a1_mean = mean(obj@meta_obs$col4a1_score, na.rm = TRUE),
    dc_mean = mean(obj@meta_obs$dc_score, na.rm = TRUE),
    tls_mean = mean(obj@meta_obs$tls_score, na.rm = TRUE),
    exclusion_mean = mean(obj@meta_obs$exclusion_score, na.rm = TRUE)
  )
  write.csv(stats_df, file.path(sample_table_dir, "salie_stats.csv"), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Done:", sample_name, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) tryCatch(analyze_salie(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
}

cat("\n=== Analysis Complete ===\n")
sink()
