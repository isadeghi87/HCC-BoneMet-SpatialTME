#!/usr/bin/env Rscript
# 00a_qc_normalization.R
# Quality control, filtering, and normalization
# HIGH-QUALITY PLOTTING VERSION

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")

OBJECTS_DIR <- "analysis/objects"
PLOTS_DIR <- "analysis/plots/qc"
LOGS_DIR <- "analysis/logs"

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGS_DIR, paste0("00a_qc_norm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== QC, Filtering, and Normalization ===\n")

SAMPLE_NAMES <- c("D1_27038840_HCCboneMet_seg", "SRout_27038840_HCC_tumour")

for (sample_name in SAMPLE_NAMES) {
  cat("\n--- Processing:", sample_name, "---\n")
  
  raw_file <- file.path(OBJECTS_DIR, paste0(sample_name, "_raw.rds"))
  norm_file <- file.path(OBJECTS_DIR, paste0(sample_name, "_norm_lognorm.rds"))
  
  if (!file.exists(raw_file)) {
    cat("  Raw object not found, skipping.\n")
    next
  }
  
  spata_obj <- readRDS(raw_file)
  cat("  Loaded:", nObs(spata_obj), "spots,", nGenes(spata_obj), "genes\n")
  
  # Calculate QC metrics
  cat("  Calculating QC metrics...\n")
  
  counts <- getCountMatrix(spata_obj)
  meta <- getMetaDf(spata_obj)
  
  # Total counts and genes per spot
  meta$nCount <- Matrix::colSums(counts)

meta$nFeature <- Matrix::colSums(counts > 0)
  
  # Mitochondrial percentage
  mt_genes <- grep("^MT-", rownames(counts), value = TRUE)
  if (length(mt_genes) > 0) {
    meta$percent_mt <- Matrix::colSums(counts[mt_genes, , drop = FALSE]) / meta$nCount * 100
  } else {
    meta$percent_mt <- 0
  }
  
  # Add to object
  spata_obj <- addMetaData(spata_obj, meta_df = meta[, c("barcodes", "nCount", "nFeature", "percent_mt")])
  
  # QC plots
  cat("  Generating QC plots...\n")
  
  qc_data <- data.frame(
    nCount = meta$nCount,
    nFeature = meta$nFeature,
    percent_mt = meta$percent_mt
  )
  
  p1 <- ggplot(qc_data, aes(x = nCount)) +
    geom_histogram(bins = 50, fill = DISCRETE_COLORS[1], color = "white") +
    scale_x_log10() +
    HQ_THEME() +
    labs(title = paste(sample_name, "- Total Counts"), x = "Total Counts (log10)", y = "Frequency")
  
  p2 <- ggplot(qc_data, aes(x = nFeature)) +
    geom_histogram(bins = 50, fill = DISCRETE_COLORS[2], color = "white") +
    HQ_THEME() +
    labs(title = "Genes per Spot", x = "Number of Genes", y = "Frequency")
  
  p3 <- ggplot(qc_data, aes(x = percent_mt)) +
    geom_histogram(bins = 50, fill = DISCRETE_COLORS[3], color = "white") +
    HQ_THEME() +
    labs(title = "Mitochondrial %", x = "% MT", y = "Frequency")
  
  p4 <- ggplot(qc_data, aes(x = nCount, y = nFeature, color = percent_mt)) +
    geom_point(alpha = 0.3, size = 0.5) +
    scale_x_log10() +
    scale_color_gradientn(colors = GRADIENT_VIRIDIS) +
    HQ_THEME() +
    labs(title = "Counts vs Genes", x = "Total Counts", y = "Genes", color = "% MT")
  
  p_combined <- (p1 | p2) / (p3 | p4)
  save_plot(p_combined, file.path(PLOTS_DIR, paste0(sample_name, "_qc_metrics.png")), width = 14, height = 10)
  
  # Filtering
  cat("  Filtering...\n")
  n_before <- nObs(spata_obj)
  
  # Filter thresholds
  min_counts <- 100
  max_mt <- 25
  min_genes <- 50
  
  keep_spots <- meta$barcodes[meta$nCount >= min_counts & 
                               meta$percent_mt <= max_mt & 
                               meta$nFeature >= min_genes]
  
  spata_obj <- subsetByBarcodes(spata_obj, barcodes = keep_spots)
  
  n_after <- nObs(spata_obj)
  cat("    Filtered:", n_before, "->", n_after, "spots (", 
      round((n_before - n_after) / n_before * 100, 1), "% removed)\n")
  
  # Normalization
  cat("  Normalizing (log-normalization)...\n")
  spata_obj <- normalizeCountMatrix(spata_obj, method = "LogNormalize")
  
  # Save
  saveRDS(spata_obj, norm_file)
  cat("  Saved:", norm_file, "\n")
  
  rm(spata_obj, counts, meta)
  gc()
}

cat("\n=== QC and Normalization Complete ===\n")
sink()
