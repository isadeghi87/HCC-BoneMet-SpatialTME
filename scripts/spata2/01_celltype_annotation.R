#!/usr/bin/env Rscript
# 01_celltype_annotation.R
# Cell type annotation using marker-based scoring - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

PLOT_DIR <- "analysis/plots/celltype"
TABLE_DIR <- "analysis/tables/celltype"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("01_celltype_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Cell Type Annotation ===\n")

CELL_TYPE_MARKERS <- list(
  Hepatocyte = c("ALB", "APOB", "SERPINA1", "TTR", "AFP", "GPC3"),
  HCC_tumour = c("AFP", "GPC3", "EPCAM", "KRT19", "KRT7"),
  T_cell = c("CD3D", "CD3E", "CD3G"),
  CD8_T = c("CD8A", "CD8B"),
  CD4_T = c("CD4"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  NK = c("NCR1", "NCAM1", "NKG7", "KLRD1", "GNLY"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD19"),
  Plasma = c("JCHAIN", "MZB1", "SDC1", "IGHG1"),
  Monocyte = c("CD14", "LYZ", "S100A8", "S100A9"),
  Macrophage = c("CD68", "CD163", "MSR1", "MARCO"),
  M1_macro = c("IL1B", "NOS2", "TNF", "CD80"),
  M2_macro = c("MRC1", "CD163", "ARG1", "CD206"),
  DC = c("CD1C", "CLEC9A", "XCR1", "LAMP3", "ITGAX"),
  pDC = c("LILRA4", "IL3RA", "CLEC4C"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  CAF = c("ACTA2", "FAP", "POSTN", "PDPN"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "ERG"),
  LSEC = c("CLEC4G", "CLEC4M", "STAB2", "FCN3"),
  Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "SOX9"),
  Stellate = c("LRAT", "RGS5", "ACTA2", "DES")
)

calc_celltype_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

annotate_celltypes <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  
  # Get standardized output prefix
  output_prefix <- get_output_prefix(sample_name)
  tissue_label <- get_tissue_label(get_tissue_from_visium(sample_name))
  cat("Output prefix:", output_prefix, "\n")
  cat("Tissue:", tissue_label, "\n")
  
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  cat("Total genes:", length(genes), "\n")
  cat("Total spots:", ncol(expr_mtr), "\n")
  
  scores <- list()
  marker_counts <- list()
  
  for (ct in names(CELL_TYPE_MARKERS)) {
    markers <- CELL_TYPE_MARKERS[[ct]]
    found_markers <- markers[markers %in% genes]
    marker_counts[[ct]] <- length(found_markers)
    if (length(found_markers) > 0) {
      scores[[ct]] <- calc_celltype_score(expr_mtr, markers, genes)
      cat(ct, ":", length(found_markers), "/", length(markers), "markers\n")
    }
  }
  
  score_df <- as.data.frame(do.call(cbind, scores))
  rownames(score_df) <- colnames(expr_mtr)
  
  for (ct in names(scores)) {
    col_name <- paste0(tolower(gsub(" ", "_", ct)), "_score")
    obj@meta_obs[[col_name]] <- scores[[ct]][obj@meta_obs$barcodes]
  }
  
  score_matrix <- as.matrix(score_df)
  score_matrix[is.na(score_matrix)] <- 0
  max_scores <- apply(score_matrix, 1, max)
  dominant_type <- apply(score_matrix, 1, function(x) names(scores)[which.max(x)])
  dominant_type[max_scores < 0.3] <- "Unassigned"
  
  obj@meta_obs$celltype <- dominant_type[obj@meta_obs$barcodes]
  obj@meta_obs$celltype_confidence <- max_scores[obj@meta_obs$barcodes]
  
  cat("\nCell type distribution:\n")
  print(table(obj@meta_obs$celltype))
  
  immune_types <- c("T_cell", "CD8_T", "CD4_T", "Treg", "NK", "B_cell", "Plasma", 
                    "Monocyte", "Macrophage", "M1_macro", "M2_macro", "DC", "pDC")
  tumour_types <- c("Hepatocyte", "HCC_tumour")
  stromal_types <- c("Fibroblast", "CAF", "Endothelial", "LSEC", "Stellate")
  
  obj@meta_obs$celltype_broad <- case_when(
    obj@meta_obs$celltype %in% immune_types ~ "Immune",
    obj@meta_obs$celltype %in% tumour_types ~ "Tumour",
    obj@meta_obs$celltype %in% stromal_types ~ "Stromal",
    obj@meta_obs$celltype == "Cholangiocyte" ~ "Epithelial",
    TRUE ~ "Other"
  )
  
  cat("\nBroad categories:\n")
  print(table(obj@meta_obs$celltype_broad))
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  # Use standardized output prefix for folder names
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # HIGH-QUALITY Cell type spatial map
  p_celltype <- ggplot(plot_df, aes(x = x, y = y, color = celltype)) +
    geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
    scale_color_manual(values = CELLTYPE_COLORS, na.value = "#d9d9d9") +
    coord_fixed() + theme_publication() +
    labs(title = paste(tissue_label, "- Cell Types"), color = "Cell Type") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)))
  save_plot(p_celltype, file.path(sample_plot_dir, paste0(output_prefix, "_celltype_spatial.png")), width = 14, height = 10)
  
  # HIGH-QUALITY Broad category map
  p_broad <- ggplot(plot_df, aes(x = x, y = y, color = celltype_broad)) +
    geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
    scale_color_manual(values = BROAD_COLORS) +
    coord_fixed() + theme_publication() +
    labs(title = paste(tissue_label, "- Broad Categories"), color = "Category") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  save_plot(p_broad, file.path(sample_plot_dir, paste0(output_prefix, "_celltype_broad.png")), width = 12, height = 10)
  
  # Key cell type score maps
  key_scores <- c("hcc_tumour_score", "cd8_t_score", "macrophage_score", "caf_score", "b_cell_score")
  key_titles <- c("HCC Tumour", "CD8+ T cells", "Macrophages", "CAF", "B cells")
  
  plots <- list()
  for (i in seq_along(key_scores)) {
    if (key_scores[i] %in% names(plot_df) && !all(is.na(plot_df[[key_scores[i]]]))) {
      plots[[i]] <- ggplot(plot_df, aes(x = x, y = y, color = .data[[key_scores[i]]])) +
        geom_point(size = 0.2, alpha = 0.85, shape = 16, stroke = 0) +
        scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9") +
        coord_fixed() + theme_publication() +
        theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm")) +
        labs(title = key_titles[i], color = "Score")
    }
  }
  
  if (length(plots) > 0) {
    combined <- wrap_plots(plots[!sapply(plots, is.null)], ncol = 3) +
      plot_annotation(title = paste(tissue_label, "- Key Cell Type Scores"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, paste0(output_prefix, "_celltype_scores.png")), width = 18, height = 12)
  }
  
  # Use standardized output prefix for table folder
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  ct_counts <- as.data.frame(table(obj@meta_obs$celltype))
  colnames(ct_counts) <- c("celltype", "count")
  ct_counts$proportion <- ct_counts$count / sum(ct_counts$count)
  ct_counts$tissue <- tissue_label
  ct_counts$platform <- "ST"
  write.csv(ct_counts, file.path(sample_table_dir, paste0(output_prefix, "_ST_celltype_counts.csv")), row.names = FALSE)
  
  broad_counts <- as.data.frame(table(obj@meta_obs$celltype_broad))
  colnames(broad_counts) <- c("category", "count")
  broad_counts$proportion <- broad_counts$count / sum(broad_counts$count)
  broad_counts$tissue <- tissue_label
  broad_counts$platform <- "ST"
  write.csv(broad_counts, file.path(sample_table_dir, paste0(output_prefix, "_ST_celltype_broad.csv")), row.names = FALSE)
  
  marker_df <- data.frame(
    celltype = names(marker_counts),
    markers_found = unlist(marker_counts),
    markers_total = sapply(CELL_TYPE_MARKERS, length)
  )
  write.csv(marker_df, file.path(sample_table_dir, paste0(output_prefix, "_ST_marker_detection.csv")), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("\nSaved annotated object:", obj_path, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) {
    tryCatch(annotate_celltypes(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
  }
}

cat("\n=== Cell Type Annotation Complete ===\n")
sink()
