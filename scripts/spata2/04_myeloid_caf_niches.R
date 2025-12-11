#!/usr/bin/env Rscript
# 04_myeloid_caf_niches.R
# Myeloid and CAF niche analysis - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")

PLOT_DIR <- "analysis/plots/myeloid_caf"
TABLE_DIR <- "analysis/tables/myeloid_caf"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("04_myeloid_caf_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Myeloid and CAF Niche Analysis ===\n")

M1_MARKERS <- c("IL1B", "NOS2", "TNF", "CD80", "CD86")
M2_MARKERS <- c("MRC1", "MSR1", "CD163", "CD206", "ARG1")
CAF_MARKERS <- c("ACTA2", "POSTN", "FAP", "COL1A1", "COL1A2")

calc_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

analyze_myeloid_caf <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  cat("M1 markers:", paste(M1_MARKERS[M1_MARKERS %in% genes], collapse = ","), "\n")
  cat("M2 markers:", paste(M2_MARKERS[M2_MARKERS %in% genes], collapse = ","), "\n")
  cat("CAF markers:", paste(CAF_MARKERS[CAF_MARKERS %in% genes], collapse = ","), "\n")
  
  obj@meta_obs$m1_score <- calc_score(expr_mtr, M1_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$m2_score <- calc_score(expr_mtr, M2_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$caf_score <- calc_score(expr_mtr, CAF_MARKERS, genes)[obj@meta_obs$barcodes]
  
  if (!all(is.na(obj@meta_obs$m1_score)) && !all(is.na(obj@meta_obs$m2_score))) {
    obj@meta_obs$m1_m2_ratio <- log2((obj@meta_obs$m1_score + 0.1) / (obj@meta_obs$m2_score + 0.1))
  }
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  scores <- c("m1_score", "m2_score", "caf_score")
  titles <- c("M1 Macrophages", "M2 Macrophages", "CAF Score")
  
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
      plot_annotation(title = paste(sample_name, "- Myeloid & CAF"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, "myeloid_caf_combined.png"), width = 18, height = 7)
  }
  
  # M1/M2 ratio map
  if ("m1_m2_ratio" %in% names(plot_df)) {
    p_ratio <- ggplot(plot_df, aes(x = x, y = y, color = m1_m2_ratio)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_gradientn(colors = DIVERGING_COLORS, na.value = "#d9d9d9", limits = c(-3, 3), oob = scales::squish) +
      coord_fixed() + theme_publication() +
      theme(legend.key.width = unit(1.5, "cm")) +
      labs(title = paste(sample_name, "- M1/M2 Ratio"),
           subtitle = "Red = M1 dominant, Blue = M2 dominant",
           color = "log2(M1/M2)")
    save_plot(p_ratio, file.path(sample_plot_dir, "m1_m2_ratio.png"), width = 12, height = 10)
  }
  
  # CD8-high neighbourhood analysis
  if ("cd8_high" %in% names(obj@meta_obs)) {
    cat("CD8-high neighbourhood analysis:\n")
    for (score in c("m1_score", "m2_score", "caf_score")) {
      if (!all(is.na(obj@meta_obs[[score]]))) {
        high_mean <- mean(obj@meta_obs[[score]][obj@meta_obs$cd8_high == "CD8-high"], na.rm = TRUE)
        low_mean <- mean(obj@meta_obs[[score]][obj@meta_obs$cd8_high == "CD8-low"], na.rm = TRUE)
        cat("  ", score, "in CD8-high:", round(high_mean, 3), "vs CD8-low:", round(low_mean, 3), "\n")
      }
    }
  }
  
  sample_table_dir <- file.path(TABLE_DIR, sample_name)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    m1_mean = mean(obj@meta_obs$m1_score, na.rm = TRUE),
    m2_mean = mean(obj@meta_obs$m2_score, na.rm = TRUE),
    caf_mean = mean(obj@meta_obs$caf_score, na.rm = TRUE)
  )
  write.csv(stats_df, file.path(sample_table_dir, "myeloid_caf_stats.csv"), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Done:", sample_name, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) tryCatch(analyze_myeloid_caf(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
}

cat("\n=== Analysis Complete ===\n")
sink()
