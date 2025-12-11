#!/usr/bin/env Rscript
# 99c_fix_imc_coords.R
# Re-process IMC files to extract coordinates properly and generate spatial plots

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(patchwork)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

IMC_DIR <- "/data/projects/HCCbone/data/IMC/27038840_Hiemens"
OBJ_DIR <- "analysis/objects/imc"
PLOT_DIR <- "analysis/plots/imc"
TABLE_DIR <- "analysis/tables/imc"
LOG_DIR <- "analysis/logs"

dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("99c_fix_imc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Fix IMC Coordinates and Generate Spatial Plots ===\n")
cat("Time:", format(Sys.time()), "\n\n")

# IMC sample mapping
IMC_SAMPLES <- list(
  list(file = "27038840_Biopsy_Spine_1_1.txt", output_prefix = "BoneMet", roi = "ROI1", tissue_label = "Bone Metastasis"),
  list(file = "27038840_Biopsy_Spine_2_2.txt", output_prefix = "BoneMet", roi = "ROI2", tissue_label = "Bone Metastasis"),
  list(file = "27038840_Biopsy_Spine_3_3.txt", output_prefix = "BoneMet", roi = "ROI3", tissue_label = "Bone Metastasis"),
  list(file = "27038840_Biopsy_Spine_4_4.txt", output_prefix = "BoneMet", roi = "ROI4", tissue_label = "Bone Metastasis"),
  list(file = "27038840_Resection_1_1.txt", output_prefix = "PrimaryTumor", roi = "ROI1", tissue_label = "Primary HCC Tumor"),
  list(file = "27038840_Resection_2_2.txt", output_prefix = "PrimaryTumor", roi = "ROI2", tissue_label = "Primary HCC Tumor"),
  list(file = "27038840_Resection_3_3.txt", output_prefix = "PrimaryTumor", roi = "ROI3", tissue_label = "Primary HCC Tumor"),
  list(file = "27038840_Resection_4_4.txt", output_prefix = "PrimaryTumor", roi = "ROI4", tissue_label = "Primary HCC Tumor"),
  list(file = "27038840_Resection_5_5.txt", output_prefix = "PrimaryTumor", roi = "ROI5", tissue_label = "Primary HCC Tumor")
)

process_imc_with_coords <- function(s) {
  cat("\n--- Processing:", s$tissue_label, "-", s$roi, "---\n")
  
  txt_path <- file.path(IMC_DIR, s$file)
  if (!file.exists(txt_path)) {
    cat("  WARNING: File not found:", txt_path, "\n")
    return(NULL)
  }
  
  # Load data
  dat <- fread(txt_path, header = TRUE)
  cat("  Loaded:", nrow(dat), "cells x", ncol(dat), "columns\n")
  
  # Extract coordinates
  x_col <- intersect(c("X", "Location_Center_X", "X_centroid", "x"), colnames(dat))[1]
  y_col <- intersect(c("Y", "Location_Center_Y", "Y_centroid", "y"), colnames(dat))[1]
  
  coords <- NULL
  if (!is.na(x_col) && !is.na(y_col)) {
    coords <- data.frame(
      cell_id = paste0(s$output_prefix, "_", s$roi, "_", seq_len(nrow(dat))),
      x = dat[[x_col]], 
      y = dat[[y_col]]
    )
    cat("  Coordinates extracted from:", x_col, ",", y_col, "\n")
    cat("  X range:", min(coords$x), "-", max(coords$x), "\n")
    cat("  Y range:", min(coords$y), "-", max(coords$y), "\n")
  } else {
    cat("  WARNING: No coordinate columns found!\n")
    return(NULL)
  }
  
  # Extract marker expression
  meta_cols <- c("CellId", "ObjectNumber", "ImageNumber", "Location_Center_X", 
                 "Location_Center_Y", "AreaShape_Area", "X_centroid", "Y_centroid",
                 "x", "y", "area", "cell_id", "image_id", "Start_push", "End_push", 
                 "Pushes_duration", "X", "Y", "Z")
  marker_cols <- setdiff(colnames(dat), meta_cols)
  marker_cols <- marker_cols[!grepl("^(Metadata|Location|AreaShape|BCKG)", marker_cols, ignore.case = TRUE)]
  
  expr_cols <- intersect(marker_cols, colnames(dat))
  if (length(expr_cols) == 0) {
    cat("  WARNING: No expression columns!\n")
    return(NULL)
  }
  
  expr_mtr <- as.matrix(dat[, ..expr_cols])
  rownames(expr_mtr) <- coords$cell_id
  
  # Clean and normalize
  expr_mtr[expr_mtr < 0] <- 0
  expr_norm <- asinh(expr_mtr / 5)
  
  cat("  Markers:", ncol(expr_norm), "\n")
  
  # Create output directories
  plot_dir <- file.path(PLOT_DIR, s$output_prefix)
  table_dir <- file.path(TABLE_DIR, s$output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save updated object
  result <- list(
    sample = paste0(s$output_prefix, "_IMC_", s$roi),
    expr_raw = expr_mtr,
    expr_norm = expr_norm,
    coords = coords,
    output_prefix = s$output_prefix,
    tissue_label = s$tissue_label,
    roi = s$roi
  )
  
  obj_path <- file.path(OBJ_DIR, paste0(s$output_prefix, "_IMC_", s$roi, "_processed.rds"))
  saveRDS(result, obj_path)
  cat("  Saved object:", obj_path, "\n")
  
  # --- Generate spatial plots ---
  cat("  Generating spatial plots...\n")
  
  # Total expression spatial
  plot_df <- data.frame(coords, total_expr = rowSums(expr_norm))
  
  p_spatial <- ggplot(plot_df, aes(x = x, y = y, color = total_expr)) +
    geom_point(size = 0.1, alpha = 0.6, shape = 16, stroke = 0) +
    scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
    coord_fixed() + theme_publication() +
    labs(title = paste(s$tissue_label, "-", s$roi, "- Total Expression (IMC)"), 
         color = "Total Expr")
  save_plot(p_spatial, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_spatial_total.png")), 
            width = 10, height = 8)
  
  # Key immune markers
  # Clean up marker names for matching (remove isotope info)
  marker_names_clean <- gsub("\\(.*\\)", "", colnames(expr_norm))
  
  key_markers <- c("CD45", "CD3", "CD8a", "CD68", "CD163", "PanCK", "ki67", "PD-L1", "CD4", "FoxP3")
  
  plots <- list()
  for (m in key_markers) {
    # Find matching column
    idx <- which(grepl(paste0("^", m), marker_names_clean, ignore.case = TRUE))
    if (length(idx) > 0) {
      marker_col <- colnames(expr_norm)[idx[1]]
      plot_df$marker_val <- expr_norm[, marker_col]
      
      p <- ggplot(plot_df, aes(x = x, y = y, color = marker_val)) +
        geom_point(size = 0.08, alpha = 0.5, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm"),
              plot.title = element_text(size = 10)) +
        labs(title = m, color = "")
      plots[[m]] <- p
    }
  }
  
  if (length(plots) >= 4) {
    combined <- wrap_plots(plots[1:min(9, length(plots))], ncol = 3) +
      plot_annotation(title = paste(s$tissue_label, "-", s$roi, "- Key Markers (IMC)"),
                      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_key_markers.png")), 
              width = 15, height = 12)
  }
  
  # Cell density heatmap
  cat("  Generating cell density map...\n")
  
  # Bin cells into grid
  x_bins <- 50
  y_bins <- 50
  plot_df$x_bin <- cut(plot_df$x, breaks = x_bins, labels = FALSE)
  plot_df$y_bin <- cut(plot_df$y, breaks = y_bins, labels = FALSE)
  
  density_df <- plot_df %>%
    group_by(x_bin, y_bin) %>%
    summarise(
      cell_count = n(),
      mean_x = mean(x),
      mean_y = mean(y),
      .groups = "drop"
    )
  
  p_density <- ggplot(density_df, aes(x = mean_x, y = mean_y, fill = cell_count)) +
    geom_tile() +
    scale_fill_gradientn(colors = GRADIENT_COLORS, na.value = "white") +
    coord_fixed() + theme_publication() +
    labs(title = paste(s$tissue_label, "-", s$roi, "- Cell Density (IMC)"),
         fill = "Cells", x = "X", y = "Y")
  save_plot(p_density, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_cell_density.png")), 
            width = 10, height = 8)
  
  cat("  Done:", s$tissue_label, "-", s$roi, "\n")
  return(TRUE)
}

# Process all samples
cat("\n========================================\n")
cat("Processing IMC samples with coordinate extraction...\n")
cat("========================================\n")

for (s in IMC_SAMPLES) {
  tryCatch(
    process_imc_with_coords(s),
    error = function(e) cat("  ERROR:", e$message, "\n")
  )
}

cat("\n\n========================================\n")
cat("IMC Coordinate Fix Complete!\n")
cat("Time:", format(Sys.time()), "\n")
cat("========================================\n")

sink()
