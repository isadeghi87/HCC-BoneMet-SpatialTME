#!/usr/bin/env Rscript
# 07_nk_neighbourhood_analysis.R
# NK cell analysis and neighbourhood distances - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")

PLOT_DIR <- "analysis/plots/nk"
TABLE_DIR <- "analysis/tables/nk"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("07_nk_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== NK Cell Neighbourhood Analysis ===\n")

NK_MARKERS <- c("NCR1", "NCAM1", "CD56", "NKG7", "KLRD1", "KLRK1", "GNLY", "PRF1")
NK_ACTIVATING <- c("KLRK1", "NCR1", "NCR2", "NCR3", "CD226")
NK_INHIBITORY <- c("KLRC1", "KIR2DL1", "KIR2DL3", "TIGIT")
CYTOTOXIC_MARKERS <- c("PRF1", "GZMA", "GZMB", "GZMK", "GNLY")

calc_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

calc_min_distance <- function(coords, query_barcodes, ref_barcodes) {
  if (length(ref_barcodes) == 0) return(rep(NA, length(query_barcodes)))
  query_coords <- coords[coords$barcodes %in% query_barcodes, ]
  ref_coords <- coords[coords$barcodes %in% ref_barcodes, ]
  distances <- sapply(seq_len(nrow(query_coords)), function(i) {
    min(sqrt((ref_coords$x - query_coords$x[i])^2 + (ref_coords$y - query_coords$y[i])^2))
  })
  names(distances) <- query_coords$barcodes
  distances
}

analyze_nk <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  cat("NK markers:", paste(NK_MARKERS[NK_MARKERS %in% genes], collapse = ","), "\n")
  cat("NK activating:", paste(NK_ACTIVATING[NK_ACTIVATING %in% genes], collapse = ","), "\n")
  cat("Cytotoxic markers:", paste(CYTOTOXIC_MARKERS[CYTOTOXIC_MARKERS %in% genes], collapse = ","), "\n")
  
  obj@meta_obs$nk_score <- calc_score(expr_mtr, NK_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$nk_activating_score <- calc_score(expr_mtr, NK_ACTIVATING, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$nk_inhibitory_score <- calc_score(expr_mtr, NK_INHIBITORY, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$cytotoxic_score <- calc_score(expr_mtr, CYTOTOXIC_MARKERS, genes)[obj@meta_obs$barcodes]
  
  if (!all(is.na(obj@meta_obs$nk_score))) {
    nk_thresh <- quantile(obj@meta_obs$nk_score, 0.85, na.rm = TRUE)
    obj@meta_obs$nk_high <- ifelse(obj@meta_obs$nk_score >= nk_thresh, "NK-high", "NK-low")
    cat("NK-high spots (>85th percentile):", sum(obj@meta_obs$nk_high == "NK-high"), "\n")
  }
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  scores <- c("nk_score", "nk_activating_score", "cytotoxic_score")
  titles <- c("NK Score", "NK Activating", "Cytotoxic Score")
  
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
    combined <- wrap_plots(plots[!sapply(plots, is.null)], ncol = 3) +
      plot_annotation(title = paste(sample_name, "- NK Analysis"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, "nk_combined.png"), width = 18, height = 7)
  }
  
  # NK-high map
  if ("nk_high" %in% names(plot_df)) {
    p_nk <- ggplot(plot_df, aes(x = x, y = y, color = nk_high)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_manual(values = c("NK-high" = "#9467bd", "NK-low" = "#d9d9d9")) +
      coord_fixed() + theme_publication() +
      labs(title = paste(sample_name, "- NK-high Regions"), color = "Status") +
      guides(color = guide_legend(override.aes = list(size = 4)))
    save_plot(p_nk, file.path(sample_plot_dir, "nk_high_map.png"), width = 12, height = 10)
  }
  
  # Neighbourhood distance analysis
  cat("\nNeighbourhood distance analysis:\n")
  nk_high_barcodes <- obj@meta_obs$barcodes[obj@meta_obs$nk_high == "NK-high"]
  
  if ("cd8_high" %in% names(obj@meta_obs)) {
    cd8_high_barcodes <- obj@meta_obs$barcodes[obj@meta_obs$cd8_high == "CD8-high"]
    if (length(cd8_high_barcodes) > 0 && length(nk_high_barcodes) > 0) {
      nk_to_cd8_dist <- calc_min_distance(coords, nk_high_barcodes, cd8_high_barcodes)
      cat("NK-high to CD8-high distance: median =", round(median(nk_to_cd8_dist, na.rm = TRUE), 1), "µm\n")
    }
  }
  
  if ("tls_hotspot" %in% names(obj@meta_obs)) {
    tls_barcodes <- obj@meta_obs$barcodes[obj@meta_obs$tls_hotspot == "TLS-high"]
    if (length(tls_barcodes) > 0 && length(nk_high_barcodes) > 0) {
      nk_to_tls_dist <- calc_min_distance(coords, nk_high_barcodes, tls_barcodes)
      cat("NK-high to TLS distance: median =", round(median(nk_to_tls_dist, na.rm = TRUE), 1), "µm\n")
    }
  }
  
  # Co-localisation: NK with CD8
  if ("cd8_high" %in% names(obj@meta_obs) && "nk_high" %in% names(obj@meta_obs)) {
    contingency <- table(obj@meta_obs$nk_high, obj@meta_obs$cd8_high)
    cat("\nNK-CD8 co-localisation:\n")
    print(contingency)
  }
  
  sample_table_dir <- file.path(TABLE_DIR, sample_name)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    nk_mean = mean(obj@meta_obs$nk_score, na.rm = TRUE),
    nk_activating_mean = mean(obj@meta_obs$nk_activating_score, na.rm = TRUE),
    cytotoxic_mean = mean(obj@meta_obs$cytotoxic_score, na.rm = TRUE),
    nk_high_count = sum(obj@meta_obs$nk_high == "NK-high", na.rm = TRUE)
  )
  write.csv(stats_df, file.path(sample_table_dir, "nk_stats.csv"), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Done:", sample_name, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) tryCatch(analyze_nk(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
}

cat("\n=== Analysis Complete ===\n")
sink()
