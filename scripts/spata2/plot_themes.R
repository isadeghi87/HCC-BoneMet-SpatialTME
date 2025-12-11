# plot_themes.R
# High-quality plotting themes and color palettes for HCCbone analysis

# High-contrast discrete color palette (20 colors)
DISCRETE_COLORS <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",

  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"
)

# Cell type colors (high contrast)
CELLTYPE_COLORS <- c(
  "Hepatocyte" = "#1f77b4", "HCC_tumour" = "#d62728",
  "T_cell" = "#2ca02c", "CD8_T" = "#17becf", "CD4_T" = "#98df8a", "Treg" = "#bcbd22",
  "NK" = "#9467bd", "B_cell" = "#ff7f0e", "Plasma" = "#ffbb78",
  "Monocyte" = "#8c564b", "Macrophage" = "#c49c94", "M1_macro" = "#e377c2", "M2_macro" = "#f7b6d2",
  "DC" = "#7f7f7f", "pDC" = "#c7c7c7",
  "Fibroblast" = "#aec7e8", "CAF" = "#ff9896", "Endothelial" = "#c5b0d5",
  "LSEC" = "#dbdb8d", "Stellate" = "#9edae5", "Cholangiocyte" = "#393b79",
  "Unassigned" = "#d9d9d9"
)

# Broad category colors
BROAD_COLORS <- c(
  "Tumour" = "#d62728", "Immune" = "#2ca02c", "Stromal" = "#1f77b4",
  "Epithelial" = "#ff7f0e", "Other" = "#d9d9d9"
)

# High-low binary colors
BINARY_COLORS <- c("high" = "#d62728", "low" = "#d9d9d9")

# Sequential gradient (for expression scores) - Turbo palette
GRADIENT_COLORS <- c("#30123b", "#4662d7", "#36aaf9", "#1ae4b6", "#72fe5e",
                      "#c8ef34", "#faba39", "#f66b19", "#ca2a04", "#7a0403")

# Diverging gradient (for ratios, correlations)
DIVERGING_COLORS <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                       "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

# Publication-quality theme
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "#f0f0f0", color = "black"),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

# Spatial plot function with high-quality settings
spatial_feature_plot <- function(df, x_col = "x", y_col = "y", color_col, 
                                  title = "", point_size = 0.3, alpha = 0.9,
                                  gradient_type = "sequential") {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = alpha, shape = 16, stroke = 0) +
    coord_fixed() +
    theme_publication() +
    labs(title = title, color = color_col)
  
  if (gradient_type == "sequential") {
    p <- p + scale_color_gradientn(colors = GRADIENT_COLORS, na.value = "#d9d9d9")
  } else if (gradient_type == "diverging") {
    p <- p + scale_color_gradientn(colors = DIVERGING_COLORS, na.value = "#d9d9d9")
  }
  
  return(p)
}

# Spatial categorical plot
spatial_category_plot <- function(df, x_col = "x", y_col = "y", color_col,
                                   title = "", point_size = 0.3, alpha = 0.9,
                                   colors = NULL) {
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = alpha, shape = 16, stroke = 0) +
    coord_fixed() +
    theme_publication() +
    labs(title = title, color = color_col)
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors, na.value = "#d9d9d9")
  } else {
    p <- p + scale_color_manual(values = DISCRETE_COLORS, na.value = "#d9d9d9")
  }
  
  return(p)
}

# High-DPI save function
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 400) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
}

cat("Plot themes and palettes loaded.\n")
