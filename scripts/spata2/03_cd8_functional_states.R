#!/usr/bin/env Rscript
# 03_cd8_functional_states.R
# CD8 functional state analysis - HIGH-QUALITY PLOTTING

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

PLOT_DIR <- "analysis/plots/cd8_functional"
TABLE_DIR <- "analysis/tables/cd8_functional"
LOG_DIR <- "analysis/logs"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("03_cd8_functional_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== CD8 Functional States Analysis ===\n")

EFFECTOR_MARKERS <- c("IFNG", "GZMB", "PRF1")
EXHAUSTION_MARKERS <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX")
TRM_MARKERS <- c("ITGAE", "CD69", "CXCR6")

calc_score <- function(expr_mtr, markers, genes) {
  m <- markers[markers %in% genes]
  if (length(m) == 0) return(rep(NA, ncol(expr_mtr)))
  colMeans(as.matrix(expr_mtr[m, , drop = FALSE]))
}

analyze_cd8_functional <- function(obj_path, sample_name) {
  cat("\n--- Processing:", sample_name, "---\n")
  output_prefix <- get_output_prefix(sample_name)
  tissue_label <- get_tissue_label(get_tissue_from_visium(sample_name))
  cat("Output prefix:", output_prefix, "\n")
  obj <- readRDS(obj_path)
  genes <- rownames(getMatrix(obj))
  expr_mtr <- getMatrix(obj, mtr_name = "LogNormalize")
  
  all_markers <- c(EFFECTOR_MARKERS, EXHAUSTION_MARKERS, TRM_MARKERS)
  available <- all_markers[all_markers %in% genes]
  missing <- all_markers[!all_markers %in% genes]
  cat("Markers available:", paste(available, collapse = ", "), "\n")
  if (length(missing) > 0) cat("Markers missing:", paste(missing, collapse = ", "), "\n")
  
  obj@meta_obs$effector_score <- calc_score(expr_mtr, EFFECTOR_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$exhaustion_score <- calc_score(expr_mtr, EXHAUSTION_MARKERS, genes)[obj@meta_obs$barcodes]
  obj@meta_obs$trm_score <- calc_score(expr_mtr, TRM_MARKERS, genes)[obj@meta_obs$barcodes]
  
  # Add individual markers
  for (marker in available) {
    obj@meta_obs[[marker]] <- as.numeric(expr_mtr[marker, obj@meta_obs$barcodes])
  }
  
  # Correlations
  if ("PDCD1" %in% available && "CD69" %in% available) {
    cor_pdcd1_cd69 <- cor(obj@meta_obs$PDCD1, obj@meta_obs$CD69, use = "complete.obs")
    cat("Correlation PDCD1 vs CD69:", round(cor_pdcd1_cd69, 3), "\n")
  }
  
  coords <- getCoordsDf(obj)
  plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
  
  sample_plot_dir <- file.path(PLOT_DIR, output_prefix)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Signature scores
  scores <- c("effector_score", "exhaustion_score", "trm_score")
  titles <- c("Effector (IFNG/GZMB/PRF1)", "Exhaustion (PD1/TIM3/LAG3)", "TRM (ITGAE/CD69/CXCR6)")
  
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
      plot_annotation(title = paste(tissue_label, "- CD8 Functional States"),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    save_plot(combined, file.path(sample_plot_dir, paste0(output_prefix, "_ST_functional_states_combined.png")), width = 18, height = 7)
  }
  
  # Effector vs Exhaustion ratio
  if (!all(is.na(plot_df$effector_score)) && !all(is.na(plot_df$exhaustion_score))) {
    plot_df$eff_exh_ratio <- log2((plot_df$effector_score + 0.1) / (plot_df$exhaustion_score + 0.1))
    
    p_ratio <- ggplot(plot_df, aes(x = x, y = y, color = eff_exh_ratio)) +
      geom_point(size = 0.3, alpha = 0.85, shape = 16, stroke = 0) +
      scale_color_gradientn(colors = DIVERGING_COLORS, na.value = "#d9d9d9", limits = c(-3, 3), oob = scales::squish) +
      coord_fixed() + theme_publication() +
      theme(legend.key.width = unit(1.5, "cm")) +
      labs(title = paste(tissue_label, "- Effector/Exhaustion Ratio"), 
           subtitle = "Red = Effector dominant, Blue = Exhaustion dominant",
           color = "log2(Eff/Exh)")
    save_plot(p_ratio, file.path(sample_plot_dir, paste0(output_prefix, "_ST_effector_exhaustion_ratio.png")), width = 12, height = 10)
  }
  
  sample_table_dir <- file.path(TABLE_DIR, output_prefix)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  stats_df <- data.frame(
    sample = sample_name,
    effector_mean = mean(obj@meta_obs$effector_score, na.rm = TRUE),
    exhaustion_mean = mean(obj@meta_obs$exhaustion_score, na.rm = TRUE),
    trm_mean = mean(obj@meta_obs$trm_score, na.rm = TRUE)
  )
  write.csv(stats_df, file.path(sample_table_dir, paste0(output_prefix, "_ST_signature_stats.csv")), row.names = FALSE)
  
  saveRDS(obj, obj_path)
  cat("Done:", sample_name, "\n")
}

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) {
    tryCatch(analyze_cd8_functional(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
  }
}

cat("\n=== Analysis Complete ===\n")
sink()
