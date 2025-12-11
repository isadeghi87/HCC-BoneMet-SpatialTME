#!/usr/bin/env Rscript
# 99b_reexport_imc_outputs.R
# Re-export IMC plots and tables with standardized naming
# Also renames existing object files to new convention

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(patchwork)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

# --- Directories ---
OBJ_DIR <- "analysis/objects/imc"
PLOT_DIR <- "analysis/plots/imc"
TABLE_DIR <- "analysis/tables/imc"
LOG_DIR <- "analysis/logs"

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("99b_reexport_imc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Re-export IMC with Standardized Naming ===\n")
cat("Time:", format(Sys.time()), "\n\n")

print_sample_mapping()

# --- IMC sample mapping ---
IMC_SAMPLES <- list(
  # Bone Metastasis ROIs
  list(old_name = "27038840_Biopsy_Spine_1_1", output_prefix = "BoneMet", roi = "ROI1", tissue_label = "Bone Metastasis"),
  list(old_name = "27038840_Biopsy_Spine_2_2", output_prefix = "BoneMet", roi = "ROI2", tissue_label = "Bone Metastasis"),
  list(old_name = "27038840_Biopsy_Spine_3_3", output_prefix = "BoneMet", roi = "ROI3", tissue_label = "Bone Metastasis"),
  list(old_name = "27038840_Biopsy_Spine_4_4", output_prefix = "BoneMet", roi = "ROI4", tissue_label = "Bone Metastasis"),
  # Primary Tumor ROIs
  list(old_name = "27038840_Resection_1_1", output_prefix = "PrimaryTumor", roi = "ROI1", tissue_label = "Primary HCC Tumor"),
  list(old_name = "27038840_Resection_2_2", output_prefix = "PrimaryTumor", roi = "ROI2", tissue_label = "Primary HCC Tumor"),
  list(old_name = "27038840_Resection_3_3", output_prefix = "PrimaryTumor", roi = "ROI3", tissue_label = "Primary HCC Tumor"),
  list(old_name = "27038840_Resection_4_4", output_prefix = "PrimaryTumor", roi = "ROI4", tissue_label = "Primary HCC Tumor"),
  list(old_name = "27038840_Resection_5_5", output_prefix = "PrimaryTumor", roi = "ROI5", tissue_label = "Primary HCC Tumor")
)

# --- Clean old directories ---
clean_old_imc_dirs <- function() {
  cat("\n--- Cleaning old IMC output directories ---\n")
  
  old_dirs <- c(
    file.path(PLOT_DIR, "27038840_Biopsy_Spine_1_1"),
    file.path(PLOT_DIR, "27038840_Biopsy_Spine_2_2"),
    file.path(PLOT_DIR, "27038840_Biopsy_Spine_3_3"),
    file.path(PLOT_DIR, "27038840_Biopsy_Spine_4_4"),
    file.path(PLOT_DIR, "27038840_Resection_1_1"),
    file.path(PLOT_DIR, "27038840_Resection_2_2"),
    file.path(PLOT_DIR, "27038840_Resection_3_3"),
    file.path(PLOT_DIR, "27038840_Resection_4_4"),
    file.path(PLOT_DIR, "27038840_Resection_5_5"),
    file.path(TABLE_DIR, "27038840_Biopsy_Spine_1_1"),
    file.path(TABLE_DIR, "27038840_Biopsy_Spine_2_2"),
    file.path(TABLE_DIR, "27038840_Biopsy_Spine_3_3"),
    file.path(TABLE_DIR, "27038840_Biopsy_Spine_4_4"),
    file.path(TABLE_DIR, "27038840_Resection_1_1"),
    file.path(TABLE_DIR, "27038840_Resection_2_2"),
    file.path(TABLE_DIR, "27038840_Resection_3_3"),
    file.path(TABLE_DIR, "27038840_Resection_4_4"),
    file.path(TABLE_DIR, "27038840_Resection_5_5")
  )
  
  for (d in old_dirs) {
    if (dir.exists(d)) {
      unlink(d, recursive = TRUE)
      cat("  Removed:", d, "\n")
    }
  }
}

# --- Rename object files ---
rename_imc_objects <- function() {
  cat("\n--- Renaming IMC object files ---\n")
  
  for (s in IMC_SAMPLES) {
    old_path <- file.path(OBJ_DIR, paste0(s$old_name, "_processed.rds"))
    new_path <- file.path(OBJ_DIR, paste0(s$output_prefix, "_IMC_", s$roi, "_processed.rds"))
    
    if (file.exists(old_path) && !file.exists(new_path)) {
      file.copy(old_path, new_path)
      cat("  Copied:", basename(old_path), "->", basename(new_path), "\n")
    } else if (file.exists(new_path)) {
      cat("  Already exists:", basename(new_path), "\n")
    }
  }
}

# --- Export IMC sample ---
export_imc_sample <- function(s) {
  cat("\n--- Processing:", s$tissue_label, "-", s$roi, "---\n")
  
  # Try new name first, then old name
  obj_path <- file.path(OBJ_DIR, paste0(s$output_prefix, "_IMC_", s$roi, "_processed.rds"))
  if (!file.exists(obj_path)) {
    obj_path <- file.path(OBJ_DIR, paste0(s$old_name, "_processed.rds"))
  }
  
  if (!file.exists(obj_path)) {
    cat("  WARNING: Object not found\n")
    return(NULL)
  }
  
  obj <- readRDS(obj_path)
  cat("  Loaded:", basename(obj_path), "\n")
  cat("  Cells:", nrow(obj$expr_norm), "\n")
  cat("  Markers:", ncol(obj$expr_norm), "\n")
  
  # Create output directories
  plot_dir <- file.path(PLOT_DIR, s$output_prefix)
  table_dir <- file.path(TABLE_DIR, s$output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Marker distribution ---
  cat("  Exporting marker distribution...\n")
  n_sample <- min(5000, nrow(obj$expr_norm))
  expr_sample <- obj$expr_norm[sample(nrow(obj$expr_norm), n_sample), ]
  expr_long <- melt(expr_sample)
  colnames(expr_long) <- c("cell", "marker", "expression")
  
  p_dist <- ggplot(expr_long, aes(x = marker, y = expression)) +
    geom_violin(fill = "#4DAF4A", alpha = 0.7, scale = "width") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          panel.grid.minor = element_blank()) +
    labs(title = paste(s$tissue_label, "-", s$roi, "- Marker Distribution (IMC)"), 
         x = "Marker", y = "Expression (arcsinh)")
  save_plot(p_dist, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_marker_dist.png")), width = 14, height = 6)
  
  # --- Marker stats table ---
  marker_stats <- data.frame(
    marker = colnames(obj$expr_norm),
    mean = round(colMeans(obj$expr_norm, na.rm = TRUE), 4),
    sd = round(apply(obj$expr_norm, 2, sd, na.rm = TRUE), 4),
    min = round(apply(obj$expr_norm, 2, min, na.rm = TRUE), 4),
    max = round(apply(obj$expr_norm, 2, max, na.rm = TRUE), 4),
    tissue = s$tissue_label,
    roi = s$roi,
    platform = "IMC"
  )
  write.csv(marker_stats, file.path(table_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_marker_stats.csv")), row.names = FALSE)
  
  # --- Spatial plot ---
  if (!is.null(obj$coords)) {
    cat("  Exporting spatial overview...\n")
    plot_df <- data.frame(obj$coords, total_expr = rowSums(obj$expr_norm))
    
    p_spatial <- ggplot(plot_df, aes(x = x, y = y, color = total_expr)) +
      geom_point(size = 0.3, alpha = 0.7, shape = 16, stroke = 0) +
      scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
      coord_fixed() + theme_publication() +
      labs(title = paste(s$tissue_label, "-", s$roi, "- Spatial Overview (IMC)"), 
           color = "Total Expression")
    save_plot(p_spatial, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_spatial.png")), width = 10, height = 8)
    
    # Key marker spatial plots
    key_markers <- c("CD45", "CD3", "CD8a", "CD68", "PanCK", "Ki67")
    available_markers <- intersect(key_markers, colnames(obj$expr_norm))
    
    if (length(available_markers) > 0) {
      plots <- list()
      for (m in available_markers) {
        plot_df[[m]] <- obj$expr_norm[, m]
        plots[[m]] <- ggplot(plot_df, aes(x = x, y = y, color = .data[[m]])) +
          geom_point(size = 0.2, alpha = 0.7, shape = 16, stroke = 0) +
          scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
          coord_fixed() + theme_publication() +
          theme(legend.position = "bottom") +
          labs(title = m, color = "Expr")
      }
      
      combined <- wrap_plots(plots, ncol = 3) +
        plot_annotation(title = paste(s$tissue_label, "-", s$roi, "- Key Markers (IMC)"),
                        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
      save_plot(combined, file.path(plot_dir, paste0(s$output_prefix, "_IMC_", s$roi, "_key_markers.png")), width = 15, height = 10)
    }
  }
  
  cat("  Done:", s$tissue_label, "-", s$roi, "\n")
  return(TRUE)
}

# --- Create combined summary ---
create_imc_summary <- function() {
  cat("\n--- Creating IMC summary across all samples ---\n")
  
  all_stats <- list()
  
  for (s in IMC_SAMPLES) {
    stats_file <- file.path(TABLE_DIR, s$output_prefix, 
                            paste0(s$output_prefix, "_IMC_", s$roi, "_marker_stats.csv"))
    if (file.exists(stats_file)) {
      df <- read.csv(stats_file)
      all_stats[[paste0(s$output_prefix, "_", s$roi)]] <- df
    }
  }
  
  if (length(all_stats) > 0) {
    combined_stats <- do.call(rbind, all_stats)
    write.csv(combined_stats, file.path(TABLE_DIR, "IMC_combined_marker_stats.csv"), row.names = FALSE)
    cat("  Saved: IMC_combined_marker_stats.csv\n")
    
    # Summary by tissue
    summary_by_tissue <- combined_stats %>%
      group_by(tissue, marker) %>%
      summarise(
        mean_expr = mean(mean, na.rm = TRUE),
        sd_expr = mean(sd, na.rm = TRUE),
        n_rois = n(),
        .groups = "drop"
      )
    write.csv(summary_by_tissue, file.path(TABLE_DIR, "IMC_summary_by_tissue.csv"), row.names = FALSE)
    cat("  Saved: IMC_summary_by_tissue.csv\n")
  }
}

# --- Main execution ---
cat("\n========================================\n")
cat("Starting IMC re-export...\n")
cat("========================================\n")

# Clean old directories
clean_old_imc_dirs()

# Rename object files
rename_imc_objects()

# Process each sample
for (s in IMC_SAMPLES) {
  tryCatch(
    export_imc_sample(s),
    error = function(e) cat("  ERROR:", e$message, "\n")
  )
}

# Create combined summary
create_imc_summary()

cat("\n\n========================================\n")
cat("IMC Re-export complete!\n")
cat("Time:", format(Sys.time()), "\n")
cat("========================================\n")

sink()
