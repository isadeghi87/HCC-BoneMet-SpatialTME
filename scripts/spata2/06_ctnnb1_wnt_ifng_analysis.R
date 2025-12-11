#!/usr/bin/env Rscript
# 06_ctnnb1_wnt_ifng_analysis.R
# CTNNB1/Wnt pathway and IFNγ/antigen presentation - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

PLOT_DIR <- "analysis/plots/wnt_ifng"
TABLE_DIR <- "analysis/tables/wnt_ifng"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("06_wnt_ifng_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== CTNNB1/Wnt and IFNγ Analysis ===\n")

WNT_TARGETS <- c("AXIN2", "MYC", "CCND1", "TCF7", "LEF1", "BIRC5", "LGR5")
CTNNB1_MARKERS <- c("CTNNB1", "APC", "GSK3B")
IFNG_RESPONSE <- c("IFNG", "IRF1", "STAT1", "JAK1", "JAK2")
HLA_CLASS_I <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2")

calc_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

analyze_wnt_ifng <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  output_prefix <- get_output_prefix(sample_name)
  tissue_label <- get_tissue_label(get_tissue_from_visium(sample_name))
  cat("Output prefix:", output_prefix, "\n")
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  cat("Wnt target markers:", paste(WNT_TARGETS[WNT_TARGETS %in% genes], collapse = ","), "\n")
  cat("CTNNB1 markers:", paste(CTNNB1_MARKERS[CTNNB1_MARKERS %in% genes], collapse = ","), "\n")
  cat("IFNg response:", paste(IFNG_RESPONSE[IFNG_RESPONSE %in% genes], collapse = ","), "\n")
  cat("HLA-I:", paste(HLA_CLASS_I[HLA_CLASS_I %in% genes], collapse = ","), "\n")
  
  obj@meta_obs$wnt_score <- calc_score(expr_mtr, WNT_TARGETS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$ctnnb1_score <- calc_score(expr_mtr, CTNNB1_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$ifng_score <- calc_score(expr_mtr, IFNG_RESPONSE, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$hla_i_score <- calc_score(expr_mtr, HLA_CLASS_I, genes)[obj@meta_obs$barcodes]
  
  if (!all(is.na(obj@meta_obs$wnt_score))) {
    wnt_thresh <- quantile(obj@meta_obs$wnt_score, 0.75, na.rm = TRUE)
    obj@meta_obs$wnt_high <- ifelse(obj@meta_obs$wnt_score >= wnt_thresh, "Wnt-high", "Wnt-low")
    cat("Wnt-high spots (>75th percentile):", sum(obj@meta_obs$wnt_high == "Wnt-high"), "\n")
  }
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  scores <- c("wnt_score", "ctnnb1_score", "ifng_score", "hla_i_score")
  titles <- c("Wnt Target Score", "CTNNB1 Score", "IFNγ Response", "HLA Class I")
  
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
      plot_annotation(title = paste(tissue_label, "- Wnt/IFNγ Analysis"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, paste0(output_prefix, "_ST_wnt_ifng_combined.png")), width = 14, height = 12)
  }
  
  # Correlation analysis: Wnt vs IFNg
  if (!all(is.na(obj@meta_obs$wnt_score)) && !all(is.na(obj@meta_obs$ifng_score))) {
    cor_wnt_ifng <- cor(obj@meta_obs$wnt_score, obj@meta_obs$ifng_score, use = "complete.obs")
    cat("Wnt-IFNg correlation:", round(cor_wnt_ifng, 3), "\n")
    
    p_scatter <- ggplot(plot_df, aes(x = wnt_score, y = ifng_score)) +
      geom_point(size = 0.3, alpha = 0.2, color = "#377eb8") +
      geom_smooth(method = "lm", se = TRUE, color = "#d62728", linewidth = 1.2) +
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank()) +
      labs(title = paste(tissue_label, "- Wnt vs IFNγ (r =", round(cor_wnt_ifng, 3), ")"),
           x = "Wnt Score", y = "IFNγ Score")
    save_plot(p_scatter, file.path(sample_plot_dir, paste0(output_prefix, "_ST_wnt_ifng_correlation.png")), width = 8, height = 7)
  }
  
  # Wnt-high zone characteristics
  if ("wnt_high" %in% names(obj@meta_obs)) {
    cat("\nWnt-high zone characteristics:\n")
    if ("cd8_score" %in% names(obj@meta_obs)) {
      wnt_high_cd8 <- mean(obj@meta_obs$cd8_score[obj@meta_obs$wnt_high == "Wnt-high"], na.rm = TRUE)
      wnt_low_cd8 <- mean(obj@meta_obs$cd8_score[obj@meta_obs$wnt_high == "Wnt-low"], na.rm = TRUE)
      cat("CD8 in Wnt-high:", round(wnt_high_cd8, 3), "vs Wnt-low:", round(wnt_low_cd8, 3), "\n")
    }
    if (!all(is.na(obj@meta_obs$ifng_score))) {
      wnt_high_ifng <- mean(obj@meta_obs$ifng_score[obj@meta_obs$wnt_high == "Wnt-high"], na.rm = TRUE)
      wnt_low_ifng <- mean(obj@meta_obs$ifng_score[obj@meta_obs$wnt_high == "Wnt-low"], na.rm = TRUE)
      cat("IFNg in Wnt-high:", round(wnt_high_ifng, 3), "vs Wnt-low:", round(wnt_low_ifng, 3), "\n")
    }
    if (!all(is.na(obj@meta_obs$hla_i_score))) {
      wnt_high_hla <- mean(obj@meta_obs$hla_i_score[obj@meta_obs$wnt_high == "Wnt-high"], na.rm = TRUE)
      wnt_low_hla <- mean(obj@meta_obs$hla_i_score[obj@meta_obs$wnt_high == "Wnt-low"], na.rm = TRUE)
      cat("HLA-I in Wnt-high:", round(wnt_high_hla, 3), "vs Wnt-low:", round(wnt_low_hla, 3), "\n")
    }
  }
  
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    wnt_mean = mean(obj@meta_obs$wnt_score, na.rm = TRUE),
    ctnnb1_mean = mean(obj@meta_obs$ctnnb1_score, na.rm = TRUE),
    ifng_mean = mean(obj@meta_obs$ifng_score, na.rm = TRUE),
    hla_i_mean = mean(obj@meta_obs$hla_i_score, na.rm = TRUE)
  )
  write.csv(stats_df, file.path(sample_table_dir, paste0(output_prefix, "_ST_wnt_ifng_stats.csv")), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Done:", sample_name, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) tryCatch(analyze_wnt_ifng(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
}

cat("\n=== Analysis Complete ===\n")
sink()
