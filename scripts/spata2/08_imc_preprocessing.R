#!/usr/bin/env Rscript
# 08_imc_preprocessing.R
# IMC data loading, normalisation, and QC - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

IMC_DIR <- "/data/projects/HCCbone/data/IMC/27038840_Hiemens"
PLOT_DIR <- "analysis/plots/imc"
TABLE_DIR <- "analysis/tables/imc"
LOG_DIR <- "analysis/logs"
OBJ_DIR <- "analysis/objects/imc"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("08_imc_preprocess_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== IMC Preprocessing ===\n")
cat("\n--- Sample Mapping ---\n")
print_sample_mapping()

txt_files <- list.files(IMC_DIR, pattern = "\\.txt$", full.names = TRUE)
cat("\nFound", length(txt_files), "TXT files:\n")
for (f in txt_files) {
  sample_base <- gsub("_\\d+_\\d+\\.txt$", "", basename(f))
  tissue <- get_tissue_from_imc(sample_base)
  tissue_label <- if (!is.na(tissue)) get_tissue_label(tissue) else "Unknown"
  cat(sprintf("  %s -> %s\n", basename(f), tissue_label))
}

process_imc_txt <- function(txt_path) {
  sample_name <- gsub("\\.txt$", "", basename(txt_path))
  sample_base <- gsub("_\\d+_\\d+$", "", sample_name)
  tissue_type <- get_tissue_from_imc(sample_base)
  tissue_label <- if (!is.na(tissue_type)) get_tissue_label(tissue_type) else "Unknown"
  output_prefix <- get_output_prefix(sample_base)
  
  cat("\n--- Processing:", sample_name, "---\n")
  cat("Tissue type:", tissue_label, "\n")
  cat("Output prefix:", output_prefix, "\n")
  cat("Matched Visium sample:", get_matched_visium(sample_base), "\n")
  
  dat <- tryCatch(fread(txt_path, header = TRUE), error = function(e) { cat("Error:", e$message, "\n"); NULL })
  if (is.null(dat)) return(NULL)
  
  cat("Dimensions:", nrow(dat), "cells x", ncol(dat), "features\n")
  cat("Columns:", paste(head(colnames(dat), 20), collapse = ", "), "...\n")
  
  meta_cols <- c("CellId", "ObjectNumber", "ImageNumber", "Location_Center_X", 
                 "Location_Center_Y", "AreaShape_Area", "X_centroid", "Y_centroid",
                 "x", "y", "area", "cell_id", "image_id", "Start_push", "End_push", 
                 "Pushes_duration", "X", "Y", "Z")
  marker_cols <- setdiff(colnames(dat), meta_cols)
  marker_cols <- marker_cols[!grepl("^(Metadata|Location|AreaShape)", marker_cols, ignore.case = TRUE)]
  
  cat("Detected marker columns:", length(marker_cols), "\n")
  
  expr_cols <- intersect(marker_cols, colnames(dat))
  if (length(expr_cols) == 0) { cat("Warning: No expression columns!\n"); return(NULL) }
  
  expr_mtr <- as.matrix(dat[, ..expr_cols])
  rownames(expr_mtr) <- paste0(sample_name, "_", seq_len(nrow(expr_mtr)))
  
  cat("\nQC metrics:\n")
  cat("  Total cells:", nrow(expr_mtr), "\n")
  cat("  Markers:", ncol(expr_mtr), "\n")
  
  n_neg <- sum(expr_mtr < 0, na.rm = TRUE)
  cat("  Negative values:", n_neg, "\n")
  expr_mtr[expr_mtr < 0] <- 0
  
  cat("  Applying arcsinh transformation (cofactor=5)...\n")
  expr_norm <- asinh(expr_mtr / 5)
  
  marker_means <- colMeans(expr_norm, na.rm = TRUE)
  marker_sds <- apply(expr_norm, 2, sd, na.rm = TRUE)
  
  stats_df <- data.frame(
    marker = colnames(expr_norm),
    mean = round(marker_means, 4),
    sd = round(marker_sds, 4),
    min = round(apply(expr_norm, 2, min, na.rm = TRUE), 4),
    max = round(apply(expr_norm, 2, max, na.rm = TRUE), 4)
  )
  
  # Extract ROI number from sample_name (e.g., "27038840_Biopsy_Spine_1_1" -> "ROI1")
  roi_num <- gsub(".*_(\\d+)_\\d+$", "\\1", sample_name)
  roi_label <- paste0("ROI", roi_num)
  
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(stats_df, file.path(sample_table_dir, paste0(output_prefix, "_IMC_", roi_label, "_marker_stats.csv")), row.names = FALSE)
  
  # Coordinates
  x_col <- intersect(c("X", "Location_Center_X", "X_centroid", "x"), colnames(dat))[1]
  y_col <- intersect(c("Y", "Location_Center_Y", "Y_centroid", "y"), colnames(dat))[1]
  
  coords <- NULL
  if (!is.na(x_col) && !is.na(y_col)) {
    coords <- data.frame(cell_id = rownames(expr_norm), x = dat[[x_col]], y = dat[[y_col]])
    cat("  Coordinates extracted from:", x_col, ",", y_col, "\n")
  }
  
  result <- list(sample = sample_name, expr_raw = expr_mtr, expr_norm = expr_norm, coords = coords, 
                 marker_stats = stats_df, output_prefix = output_prefix, tissue_label = tissue_label, roi = roi_label)
  saveRDS(result, file.path(OBJ_DIR, paste0(output_prefix, "_IMC_", roi_label, "_processed.rds")))
  cat("  Saved:", file.path(OBJ_DIR, paste0(output_prefix, "_IMC_", roi_label, "_processed.rds")), "\n")
  
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # HIGH-QUALITY marker distribution plot
  n_sample <- min(5000, nrow(expr_norm))
  expr_sample <- expr_norm[sample(nrow(expr_norm), n_sample), ]
  expr_long <- melt(expr_sample)
  colnames(expr_long) <- c("cell", "marker", "expression")
  
  p_dist <- ggplot(expr_long, aes(x = marker, y = expression)) +
    geom_violin(fill = "#4DAF4A", alpha = 0.7, scale = "width") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid.minor = element_blank()) +
    labs(title = paste(tissue_label, "-", roi_label, "- Marker Distribution (arcsinh)"), x = "Marker", y = "Expression")
  save_plot(p_dist, file.path(sample_plot_dir, paste0(output_prefix, "_IMC_", roi_label, "_marker_dist.png")), width = 14, height = 6)
  
  # Spatial plot
  if (!is.null(coords)) {
    sample_expr <- data.frame(coords, total_expr = rowSums(expr_norm))
    p_spatial <- ggplot(sample_expr, aes(x = x, y = y, color = total_expr)) +
      geom_point(size = 0.1, alpha = 0.6, shape = 16, stroke = 0) +
      scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
      coord_fixed() + theme_publication() +
      labs(title = paste(tissue_label, "-", roi_label, "- Spatial Distribution"), color = "Total Expr")
    save_plot(p_spatial, file.path(sample_plot_dir, paste0(output_prefix, "_IMC_", roi_label, "_spatial.png")), width = 10, height = 8)
  }
  
  cat("Done:", sample_name, "\n")
  return(result)
}

results <- list()
for (txt_file in txt_files) {
  res <- process_imc_txt(txt_file)
  if (!is.null(res)) results[[res$sample]] <- res
}

cat("\n=== Preprocessing Complete ===\n")
cat("Processed", length(results), "samples\n")

sink()
