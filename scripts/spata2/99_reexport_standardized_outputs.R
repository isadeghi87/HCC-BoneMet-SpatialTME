#!/usr/bin/env Rscript
# 99_reexport_standardized_outputs.R
# Re-export all plots and tables with standardized naming
# This script loads existing objects and regenerates outputs without re-running analysis

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

# --- Directories ---
OBJ_DIR <- "analysis/objects"
BASE_PLOT_DIR <- "analysis/plots"
BASE_TABLE_DIR <- "analysis/tables"
LOG_DIR <- "analysis/logs"

log_file <- file.path(LOG_DIR, paste0("99_reexport_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Re-export with Standardized Naming ===\n")
cat("Time:", format(Sys.time()), "\n\n")

print_sample_mapping()

# --- Sample definitions ---
SAMPLES <- list(
  list(
    obj_path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds",
    sample_name = "D1_27038840_HCCboneMet_seg",
    output_prefix = "BoneMet",
    tissue_label = "Bone Metastasis"
  ),
  list(
    obj_path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds",
    sample_name = "SRout_27038840_HCC_tumour",
    output_prefix = "PrimaryTumor",
    tissue_label = "Primary HCC Tumor"
  )
)

# --- Clean old directories ---
clean_old_dirs <- function() {
  cat("\n--- Cleaning old output directories ---\n")
  
  old_dirs <- c(
    # Old ST directories
    file.path(BASE_PLOT_DIR, "celltype/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "celltype/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "cd8_spatial/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "cd8_spatial/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "cd8_functional/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "cd8_functional/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "myeloid_caf/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "myeloid_caf/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "salie/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "salie/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "wnt_ifng/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "wnt_ifng/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "nk/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "nk/SRout_27038840_HCC_tumour"),
    file.path(BASE_PLOT_DIR, "spatial_interactions/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_PLOT_DIR, "spatial_interactions/SRout_27038840_HCC_tumour"),
    # Old table directories
    file.path(BASE_TABLE_DIR, "celltype/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "celltype/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "cd8_spatial/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "cd8_spatial/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "cd8_functional/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "cd8_functional/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "myeloid_caf/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "myeloid_caf/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "salie/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "salie/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "wnt_ifng/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "wnt_ifng/SRout_27038840_HCC_tumour"),
    file.path(BASE_TABLE_DIR, "nk/D1_27038840_HCCboneMet_seg"),
    file.path(BASE_TABLE_DIR, "nk/SRout_27038840_HCC_tumour"),
    # Old IMC directories
    file.path(BASE_PLOT_DIR, "imc/27038840_Biopsy_Spine_1_1"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Biopsy_Spine_2_2"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Biopsy_Spine_3_3"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Biopsy_Spine_4_4"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Resection_1_1"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Resection_2_2"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Resection_3_3"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Resection_4_4"),
    file.path(BASE_PLOT_DIR, "imc/27038840_Resection_5_5")
  )
  
  for (d in old_dirs) {
    if (dir.exists(d)) {
      unlink(d, recursive = TRUE)
      cat("  Removed:", d, "\n")
    }
  }
}

# --- Export cell type annotation ---
export_celltype <- function(obj, output_prefix, tissue_label) {
  cat("\n[1/7] Exporting cell type annotation...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "celltype", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "celltype", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  if ("celltype" %in% names(plot_df)) {
    # Cell type spatial map
    p_celltype <- ggplot(plot_df, aes(x = x, y = y, color = celltype)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_manual(values = CELLTYPE_COLORS, na.value = "#d9d9d9") +
      coord_fixed() + theme_publication() +
      labs(title = paste(tissue_label, "- Cell Types (ST)"), color = "Cell Type") +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)))
    save_plot(p_celltype, file.path(plot_dir, paste0(output_prefix, "_ST_celltype_spatial.png")), width = 14, height = 10)
    
    # Cell type counts
    ct_counts <- as.data.frame(table(obj@meta_obs$celltype))
    colnames(ct_counts) <- c("celltype", "count")
    ct_counts$proportion <- ct_counts$count / sum(ct_counts$count)
    ct_counts$tissue <- tissue_label
    ct_counts$platform <- "ST"
    write.csv(ct_counts, file.path(table_dir, paste0(output_prefix, "_ST_celltype_counts.csv")), row.names = FALSE)
    cat("  Saved celltype outputs\n")
  }
  
  if ("celltype_broad" %in% names(plot_df)) {
    p_broad <- ggplot(plot_df, aes(x = x, y = y, color = celltype_broad)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_manual(values = BROAD_COLORS) +
      coord_fixed() + theme_publication() +
      labs(title = paste(tissue_label, "- Broad Categories (ST)"), color = "Category") +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
    save_plot(p_broad, file.path(plot_dir, paste0(output_prefix, "_ST_celltype_broad.png")), width = 12, height = 10)
  }
}

# --- Export CD8 spatial ---
export_cd8_spatial <- function(obj, output_prefix, tissue_label) {
  cat("\n[2/7] Exporting CD8 spatial enrichment...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "cd8_spatial", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "cd8_spatial", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  cd8_cols <- grep("cd8|CD8", names(plot_df), value = TRUE)
  for (col in cd8_cols[1:min(3, length(cd8_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved CD8 spatial outputs\n")
}

# --- Export CD8 functional ---
export_cd8_functional <- function(obj, output_prefix, tissue_label) {
  cat("\n[3/7] Exporting CD8 functional states...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "cd8_functional", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "cd8_functional", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  func_cols <- grep("exhaust|cytotox|memory|naive|effector", names(plot_df), value = TRUE, ignore.case = TRUE)
  for (col in func_cols[1:min(4, length(func_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved CD8 functional outputs\n")
}

# --- Export myeloid/CAF ---
export_myeloid_caf <- function(obj, output_prefix, tissue_label) {
  cat("\n[4/7] Exporting myeloid/CAF niches...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "myeloid_caf", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "myeloid_caf", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  myeloid_cols <- grep("macro|m1_|m2_|caf|tam|myeloid", names(plot_df), value = TRUE, ignore.case = TRUE)
  for (col in myeloid_cols[1:min(4, length(myeloid_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved myeloid/CAF outputs\n")
}

# --- Export Salié immunotypes ---
export_salie <- function(obj, output_prefix, tissue_label) {
  cat("\n[5/7] Exporting Salié immunotypes...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "salie", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "salie", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  salie_cols <- grep("salie|immunotype|inflam|desert|exclude", names(plot_df), value = TRUE, ignore.case = TRUE)
  for (col in salie_cols[1:min(4, length(salie_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      if (is.numeric(plot_df[[col]])) {
        p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
          geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
          scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
          coord_fixed() + theme_publication() +
          labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      } else {
        p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
          geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
          scale_color_manual(values = DISCRETE_COLORS) +
          coord_fixed() + theme_publication() +
          labs(title = paste(tissue_label, "-", col, "(ST)"), color = col)
      }
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved Salié outputs\n")
}

# --- Export Wnt/IFNg ---
export_wnt_ifng <- function(obj, output_prefix, tissue_label) {
  cat("\n[6/7] Exporting Wnt/IFNg pathway...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "wnt_ifng", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "wnt_ifng", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  wnt_cols <- grep("wnt|ctnnb1|ifn|antigen|beta.catenin", names(plot_df), value = TRUE, ignore.case = TRUE)
  for (col in wnt_cols[1:min(4, length(wnt_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved Wnt/IFNg outputs\n")
}

# --- Export NK ---
export_nk <- function(obj, output_prefix, tissue_label) {
  cat("\n[7/7] Exporting NK cell analysis...\n")
  
  plot_dir <- file.path(BASE_PLOT_DIR, "nk", output_prefix)
  table_dir <- file.path(BASE_TABLE_DIR, "nk", output_prefix)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  nk_cols <- grep("nk|ncr1|ncam|nkg7|klr", names(plot_df), value = TRUE, ignore.case = TRUE)
  for (col in nk_cols[1:min(4, length(nk_cols))]) {
    if (!all(is.na(plot_df[[col]]))) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        labs(title = paste(tissue_label, "-", col, "(ST)"), color = "Score")
      save_plot(p, file.path(plot_dir, paste0(output_prefix, "_ST_", gsub("[^A-Za-z0-9]", "_", col), ".png")), width = 12, height = 10)
    }
  }
  cat("  Saved NK outputs\n")
}

# --- Main execution ---
cat("\n========================================\n")
cat("Starting re-export...\n")
cat("========================================\n")

# Clean old directories
clean_old_dirs()

# Process each sample
for (s in SAMPLES) {
  cat("\n\n========================================\n")
  cat("Processing:", s$tissue_label, "\n")
  cat("Output prefix:", s$output_prefix, "\n")
  cat("========================================\n")
  
  if (!file.exists(s$obj_path)) {
    cat("WARNING: Object not found:", s$obj_path, "\n")
    next
  }
  
  obj <- readRDS(s$obj_path)
  cat("Loaded object:", s$obj_path, "\n")
  cat("Spots:", nrow(obj@meta_obs), "\n")
  cat("Available metadata columns:", paste(head(names(obj@meta_obs), 20), collapse = ", "), "...\n")
  
  # Export all analyses
  tryCatch(export_celltype(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in celltype:", e$message, "\n"))
  tryCatch(export_cd8_spatial(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in cd8_spatial:", e$message, "\n"))
  tryCatch(export_cd8_functional(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in cd8_functional:", e$message, "\n"))
  tryCatch(export_myeloid_caf(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in myeloid_caf:", e$message, "\n"))
  tryCatch(export_salie(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in salie:", e$message, "\n"))
  tryCatch(export_wnt_ifng(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in wnt_ifng:", e$message, "\n"))
  tryCatch(export_nk(obj, s$output_prefix, s$tissue_label), 
           error = function(e) cat("Error in nk:", e$message, "\n"))
  
  rm(obj)
  gc()
}

cat("\n\n========================================\n")
cat("Re-export complete!\n")
cat("Time:", format(Sys.time()), "\n")
cat("========================================\n")

sink()
