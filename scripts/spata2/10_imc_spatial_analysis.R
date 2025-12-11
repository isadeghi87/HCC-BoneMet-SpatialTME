#!/usr/bin/env Rscript
# 10_imc_spatial_analysis.R
# Comprehensive IMC Spatial Analysis (Bodenmiller-style)
# Includes: Interaction graphs, spatial visualization, community analysis,
# cellular neighborhoods, patch detection, and interaction statistics
# HIGH-QUALITY PLOTTING VERSION

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(igraph)
  library(dbscan)
  library(patchwork)
})

source("scripts/spata2/plot_themes.R")

OBJ_DIR <- "analysis/objects/imc"
PLOT_DIR <- "analysis/plots/imc"
TABLE_DIR <- "analysis/tables/imc"
LOG_DIR <- "analysis/logs"

log_file <- file.path(LOG_DIR, paste0("10_imc_spatial_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Comprehensive IMC Spatial Analysis ===\n")
cat("Components: Interaction graphs, Spatial visualization, Community analysis,\n")
cat("            Cellular neighborhoods, Patch detection, Interaction statistics\n\n")

# ============================================================================
# 1. SPATIAL INTERACTION GRAPH CONSTRUCTION
# ============================================================================

build_spatial_graph <- function(coords, radius = 30, k_neighbors = NULL) {
  cat("  Building spatial interaction graph...\n")
  n <- nrow(coords)
  
  edges <- list()
  
  # Radius-based (sample for efficiency)
  sample_size <- min(n, 10000)
  sample_idx <- sample(n, sample_size)
  
  for (i in sample_idx) {
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    for (j in neighbors) {
      edges[[length(edges) + 1]] <- c(i, j)
    }
  }
  
  if (length(edges) == 0) {
    cat("    No edges found\n")
    return(NULL)
  }
  
  edge_matrix <- do.call(rbind, edges)
  g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)
  g <- igraph::simplify(g)
  
  cat("    Nodes:", igraph::vcount(g), "Edges:", igraph::ecount(g), "\n")
  return(g)
}

# ============================================================================
# 2. CELLULAR NEIGHBORHOOD ANALYSIS (CN)
# ============================================================================

compute_cellular_neighborhoods <- function(coords, clusters, radius = 30) {
  cat("  Computing cellular neighborhoods...\n")
  n <- nrow(coords)
  unique_clusters <- sort(unique(clusters))
  n_clusters <- length(unique_clusters)
  
  cn_matrix <- matrix(0, nrow = n, ncol = n_clusters,
                      dimnames = list(NULL, paste0("C", unique_clusters)))
  
  sample_idx <- if (n > 15000) sample(n, 15000) else seq_len(n)
  
  for (idx in seq_along(sample_idx)) {
    i <- sample_idx[idx]
    if (idx %% 3000 == 0) cat("    Processed", idx, "/", length(sample_idx), "\n")
    
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    
    if (length(neighbors) > 0) {
      neighbor_clusters <- clusters[neighbors]
      for (cl in unique_clusters) {
        cn_matrix[i, paste0("C", cl)] <- sum(neighbor_clusters == cl)
      }
    }
  }
  
  row_totals <- rowSums(cn_matrix)
  row_totals[row_totals == 0] <- 1
  cn_matrix_norm <- cn_matrix / row_totals
  
  return(list(raw = cn_matrix, normalized = cn_matrix_norm, sample_idx = sample_idx))
}

# ============================================================================
# 3. SPATIAL COMMUNITY DETECTION
# ============================================================================

detect_spatial_communities <- function(coords, clusters, eps = 50, minPts = 5) {
  cat("  Detecting spatial communities (DBSCAN)...\n")
  
  db_result <- dbscan::dbscan(as.matrix(coords[, c("x", "y")]), eps = eps, minPts = minPts)
  
  communities <- db_result$cluster
  n_communities <- max(communities)
  
  cat("    Found", n_communities, "spatial communities\n")
  cat("    Noise points:", sum(communities == 0), "\n")
  
  community_composition <- data.frame()
  for (comm in 1:n_communities) {
    comm_cells <- which(communities == comm)
    if (length(comm_cells) > 0) {
      comm_clusters <- clusters[comm_cells]
      comp <- table(comm_clusters) / length(comm_cells)
      comp_df <- data.frame(community = comm, 
                            cluster = names(comp), 
                            proportion = as.numeric(comp))
      community_composition <- rbind(community_composition, comp_df)
    }
  }
  
  return(list(membership = communities, n_communities = n_communities,
              composition = community_composition))
}

# ============================================================================
# 4. PATCH DETECTION (Contiguous regions of same cell type)
# ============================================================================

detect_patches <- function(coords, clusters, eps = 20, minPts = 3) {
  cat("  Detecting cell type patches...\n")
  
  unique_clusters <- sort(unique(clusters))
  patch_results <- list()
  
  for (cl in unique_clusters) {
    cl_idx <- which(clusters == cl)
    if (length(cl_idx) < minPts) next
    
    cl_coords <- coords[cl_idx, c("x", "y")]
    
    db <- tryCatch({
      dbscan::dbscan(as.matrix(cl_coords), eps = eps, minPts = minPts)
    }, error = function(e) NULL)
    
    if (!is.null(db)) {
      n_patches <- max(db$cluster)
      patch_sizes <- table(db$cluster[db$cluster > 0])
      
      patch_results[[paste0("C", cl)]] <- list(
        n_patches = n_patches,
        sizes = as.numeric(patch_sizes),
        mean_size = ifelse(length(patch_sizes) > 0, mean(as.numeric(patch_sizes)), 0),
        max_size = ifelse(length(patch_sizes) > 0, max(as.numeric(patch_sizes)), 0)
      )
    }
  }
  
  cat("    Patches detected for", length(patch_results), "cell types\n")
  return(patch_results)
}

# ============================================================================
# 5. INTERACTION ANALYSIS (Enrichment/Avoidance)
# ============================================================================

compute_interaction_enrichment <- function(coords, clusters, radius = 30, n_permutations = 50) {
  cat("  Computing interaction enrichment (permutation test)...\n")
  
  unique_clusters <- sort(unique(clusters))
  n_clusters <- length(unique_clusters)
  n <- nrow(coords)
  
  observed <- matrix(0, nrow = n_clusters, ncol = n_clusters,
                     dimnames = list(paste0("C", unique_clusters), paste0("C", unique_clusters)))
  
  sample_idx <- if (n > 10000) sample(n, 10000) else seq_len(n)
  
  for (i in sample_idx) {
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    ci <- clusters[i]
    for (j in neighbors) {
      cj <- clusters[j]
      observed[paste0("C", ci), paste0("C", cj)] <- observed[paste0("C", ci), paste0("C", cj)] + 1
    }
  }
  
  cat("    Running", n_permutations, "permutations...\n")
  perm_results <- array(0, dim = c(n_clusters, n_clusters, n_permutations))
  
  for (p in 1:n_permutations) {
    if (p %% 10 == 0) cat("      Permutation", p, "\n")
    
    perm_clusters <- sample(clusters)
    perm_mat <- matrix(0, nrow = n_clusters, ncol = n_clusters)
    
    for (i in sample_idx[1:min(5000, length(sample_idx))]) {
      dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
      neighbors <- which(dists > 0 & dists <= radius)
      ci <- perm_clusters[i]
      for (j in neighbors) {
        cj <- perm_clusters[j]
        ci_idx <- which(unique_clusters == ci)
        cj_idx <- which(unique_clusters == cj)
        perm_mat[ci_idx, cj_idx] <- perm_mat[ci_idx, cj_idx] + 1
      }
    }
    perm_results[,,p] <- perm_mat
  }
  
  perm_means <- apply(perm_results, c(1,2), mean)
  perm_sds <- apply(perm_results, c(1,2), sd)
  perm_sds[perm_sds == 0] <- 1
  
  z_scores <- (observed - perm_means) / perm_sds
  p_values <- 2 * pnorm(-abs(z_scores))
  
  dimnames(z_scores) <- list(paste0("C", unique_clusters), paste0("C", unique_clusters))
  dimnames(p_values) <- list(paste0("C", unique_clusters), paste0("C", unique_clusters))
  
  return(list(observed = observed, expected = perm_means, 
              z_scores = z_scores, p_values = p_values))
}

# ============================================================================
# 6. SPATIAL CONTEXT ANALYSIS
# ============================================================================

analyze_spatial_context <- function(coords, clusters, grid_size = 50) {
  cat("  Analyzing spatial context...\n")
  
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  
  x_breaks <- seq(x_range[1], x_range[2], length.out = grid_size + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = grid_size + 1)
  
  coords$x_bin <- cut(coords$x, breaks = x_breaks, labels = FALSE, include.lowest = TRUE)
  coords$y_bin <- cut(coords$y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  
  density_grid <- coords %>%
    group_by(x_bin, y_bin) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  cluster_grid <- coords %>%
    mutate(cluster = clusters) %>%
    group_by(x_bin, y_bin, cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(x_bin, y_bin) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  density_threshold <- quantile(density_grid$n_cells, 0.9, na.rm = TRUE)
  hotspots <- density_grid %>% filter(n_cells >= density_threshold)
  
  cat("    Grid cells:", nrow(density_grid), "\n")
  cat("    Hotspots (top 10%):", nrow(hotspots), "\n")
  
  return(list(density = density_grid, cluster_composition = cluster_grid, 
              hotspots = hotspots, grid_size = grid_size))
}

# ============================================================================
# MAIN PROCESSING FUNCTION
# ============================================================================

process_imc_spatial_comprehensive <- function(sample_path) {
  sample_name <- gsub("_processed\\.rds$", "", basename(sample_path))
  cat("\n", strrep("=", 60), "\n")
  cat("Processing:", sample_name, "\n")
  cat(strrep("=", 60), "\n")
  
  dat <- readRDS(sample_path)
  
  if (is.null(dat$coords)) { cat("No coordinates, skipping\n"); return(NULL) }
  if (is.null(dat$clusters_sub)) { cat("No clusters, skipping\n"); return(NULL) }
  
  coords <- dat$coords[dat$subsample_idx, ]
  clusters <- dat$clusters_sub
  
  cat("Cells:", nrow(coords), ", Clusters:", length(unique(clusters)), "\n\n")
  
  sample_plot_dir <- file.path(PLOT_DIR, sample_name, "spatial")
  sample_table_dir <- file.path(TABLE_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  # 1. Spatial Interaction Graph
  cat("\n[1/6] Spatial Interaction Graph\n")
  tryCatch({
    spatial_graph <- build_spatial_graph(coords, radius = 30)
    results$graph <- spatial_graph
    
    if (!is.null(spatial_graph)) {
      graph_stats <- data.frame(
        metric = c("nodes", "edges", "density", "avg_degree", "components"),
        value = c(
          igraph::vcount(spatial_graph),
          igraph::ecount(spatial_graph),
          igraph::edge_density(spatial_graph),
          mean(igraph::degree(spatial_graph)),
          igraph::components(spatial_graph)$no
        )
      )
      write.csv(graph_stats, file.path(sample_table_dir, "graph_statistics.csv"), row.names = FALSE)
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 2. Cellular Neighborhood Analysis
  cat("\n[2/6] Cellular Neighborhood Analysis\n")
  tryCatch({
    cn_result <- compute_cellular_neighborhoods(coords, clusters, radius = 30)
    results$neighborhoods <- cn_result
    
    cn_df <- as.data.frame(cn_result$normalized[cn_result$sample_idx, ])
    cn_df$cell_cluster <- clusters[cn_result$sample_idx]
    write.csv(cn_df, file.path(sample_table_dir, "cellular_neighborhoods.csv"), row.names = FALSE)
    
    cn_summary <- cn_df %>%
      group_by(cell_cluster) %>%
      summarise(across(starts_with("C"), mean, na.rm = TRUE)) %>%
      pivot_longer(-cell_cluster, names_to = "neighbor_type", values_to = "proportion")
    
    p_cn <- ggplot(cn_summary, aes(x = factor(cell_cluster), y = neighbor_type, fill = proportion)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradientn(colors = GRADIENT_VIRIDIS) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(sample_name, "- Cellular Neighborhood Composition"),
           x = "Cell Type", y = "Neighbor Type", fill = "Proportion")
    save_plot(p_cn, file.path(sample_plot_dir, "cellular_neighborhoods.png"), width = 10, height = 8)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 3. Spatial Community Detection
  cat("\n[3/6] Spatial Community Detection\n")
  tryCatch({
    communities <- detect_spatial_communities(coords, clusters, eps = 50, minPts = 5)
    results$communities <- communities
    
    comm_df <- data.frame(x = coords$x, y = coords$y, cluster = clusters, community = communities$membership)
    write.csv(comm_df, file.path(sample_table_dir, "spatial_communities.csv"), row.names = FALSE)
    
    if (nrow(communities$composition) > 0) {
      write.csv(communities$composition, file.path(sample_table_dir, "community_composition.csv"), row.names = FALSE)
    }
    
    plot_df <- comm_df %>% filter(community > 0)
    if (nrow(plot_df) > 0) {
      p_comm <- ggplot(plot_df, aes(x = x, y = y, color = factor(community))) +
        geom_point(size = 0.3, alpha = 0.6) +
        scale_color_manual(values = DISCRETE_COLORS) +
        HQ_THEME() +
        coord_fixed() +
        labs(title = paste(sample_name, "- Spatial Communities"), x = "X (µm)", y = "Y (µm)", color = "Community")
      save_plot(p_comm, file.path(sample_plot_dir, "spatial_communities.png"), width = 12, height = 10)
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 4. Patch Detection
  cat("\n[4/6] Patch Detection\n")
  tryCatch({
    patches <- detect_patches(coords, clusters, eps = 20, minPts = 3)
    results$patches <- patches
    
    patch_summary <- do.call(rbind, lapply(names(patches), function(cl) {
      data.frame(cluster = cl, n_patches = patches[[cl]]$n_patches,
                 mean_size = patches[[cl]]$mean_size, max_size = patches[[cl]]$max_size)
    }))
    write.csv(patch_summary, file.path(sample_table_dir, "patch_summary.csv"), row.names = FALSE)
    
    p_patches <- ggplot(patch_summary, aes(x = reorder(cluster, -n_patches), y = n_patches, fill = mean_size)) +
      geom_col(width = 0.7) +
      scale_fill_gradientn(colors = GRADIENT_INFERNO) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(sample_name, "- Cell Type Patches"), x = "Cell Type", y = "Number of Patches", fill = "Mean Size")
    save_plot(p_patches, file.path(sample_plot_dir, "patch_detection.png"), width = 10, height = 6)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 5. Interaction Enrichment Analysis
  cat("\n[5/6] Interaction Enrichment Analysis\n")
  tryCatch({
    interactions <- compute_interaction_enrichment(coords, clusters, radius = 30, n_permutations = 50)
    results$interactions <- interactions
    
    write.csv(interactions$observed, file.path(sample_table_dir, "interaction_observed.csv"))
    write.csv(interactions$z_scores, file.path(sample_table_dir, "interaction_zscore.csv"))
    write.csv(interactions$p_values, file.path(sample_table_dir, "interaction_pvalues.csv"))
    
    z_melt <- melt(interactions$z_scores)
    colnames(z_melt) <- c("Source", "Target", "Z_score")
    z_melt$Z_score <- pmax(pmin(z_melt$Z_score, 5), -5)
    
    p_zscore <- ggplot(z_melt, aes(x = Target, y = Source, fill = Z_score)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradientn(colors = DIVERGENT_RDBU, limits = c(-5, 5)) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(sample_name, "- Interaction Enrichment"),
           subtitle = "Z-score: >0 enriched (red), <0 avoided (blue)",
           x = "Target", y = "Source", fill = "Z-score")
    save_plot(p_zscore, file.path(sample_plot_dir, "interaction_enrichment.png"), width = 10, height = 8)
    
    sig_interactions <- z_melt %>%
      left_join(melt(interactions$p_values) %>% rename(Source = Var1, Target = Var2, p_value = value), by = c("Source", "Target")) %>%
      filter(p_value < 0.05) %>%
      arrange(desc(abs(Z_score)))
    write.csv(sig_interactions, file.path(sample_table_dir, "significant_interactions.csv"), row.names = FALSE)
    cat("    Significant interactions:", nrow(sig_interactions), "\n")
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 6. Spatial Context Analysis
  cat("\n[6/6] Spatial Context Analysis\n")
  tryCatch({
    spatial_context <- analyze_spatial_context(coords, clusters, grid_size = 50)
    results$spatial_context <- spatial_context
    
    write.csv(spatial_context$density, file.path(sample_table_dir, "spatial_density.csv"), row.names = FALSE)
    write.csv(spatial_context$hotspots, file.path(sample_table_dir, "spatial_hotspots.csv"), row.names = FALSE)
    
    p_density <- ggplot(spatial_context$density, aes(x = x_bin, y = y_bin, fill = n_cells)) +
      geom_tile() +
      scale_fill_gradientn(colors = GRADIENT_INFERNO) +
      HQ_THEME() +
      coord_fixed() +
      labs(title = paste(sample_name, "- Cell Density Map"), x = "X grid", y = "Y grid", fill = "Cells")
    save_plot(p_density, file.path(sample_plot_dir, "cell_density_map.png"), width = 10, height = 8)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 7. Spatial Visualization (Cell type map)
  cat("\n[Bonus] Spatial Visualization\n")
  tryCatch({
    plot_df <- data.frame(x = coords$x, y = coords$y, cluster = factor(clusters))
    if (nrow(plot_df) > 50000) plot_df <- plot_df[sample(nrow(plot_df), 50000), ]
    
    p_spatial <- ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.2, alpha = 0.5) +
      scale_color_manual(values = DISCRETE_COLORS[1:length(unique(clusters))]) +
      HQ_THEME() +
      coord_fixed() +
      labs(title = paste(sample_name, "- Cell Type Spatial Map"), x = "X (µm)", y = "Y (µm)", color = "Cluster") +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    save_plot(p_spatial, file.path(sample_plot_dir, "celltype_spatial_map.png"), width = 14, height = 12)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  dat$spatial_analysis <- results
  saveRDS(dat, sample_path)
  
  cat("\n✓ Spatial analysis complete for", sample_name, "\n")
  return(results)
}

# ============================================================================
# PROCESS ALL SAMPLES
# ============================================================================

processed_files <- list.files(OBJ_DIR, pattern = "_processed\\.rds$", full.names = TRUE)
cat("Found", length(processed_files), "processed IMC samples\n")

all_results <- list()
for (f in processed_files) {
  result <- tryCatch(process_imc_spatial_comprehensive(f), error = function(e) { cat("FATAL ERROR:", e$message, "\n"); NULL })
  sample_name <- gsub("_processed\\.rds$", "", basename(f))
  all_results[[sample_name]] <- result
}

# Cross-sample summary
cat("\n", strrep("=", 60), "\n")
cat("Cross-Sample Summary\n")
cat(strrep("=", 60), "\n")

all_patches <- do.call(rbind, lapply(names(all_results), function(s) {
  if (!is.null(all_results[[s]]$patches)) {
    patch_df <- do.call(rbind, lapply(names(all_results[[s]]$patches), function(cl) {
      data.frame(sample = s, cluster = cl, n_patches = all_results[[s]]$patches[[cl]]$n_patches, mean_size = all_results[[s]]$patches[[cl]]$mean_size)
    }))
    return(patch_df)
  }
  return(NULL)
}))

if (!is.null(all_patches) && nrow(all_patches) > 0) {
  write.csv(all_patches, file.path(TABLE_DIR, "all_samples_patch_summary.csv"), row.names = FALSE)
  
  p_all_patches <- ggplot(all_patches, aes(x = cluster, y = n_patches, fill = sample)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = DISCRETE_COLORS) +
    HQ_THEME() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Patch Detection Across Samples", x = "Cell Type", y = "Number of Patches")
  save_plot(p_all_patches, file.path(PLOT_DIR, "all_samples_patches.png"), width = 14, height = 8)
}

cat("\n=== Comprehensive IMC Spatial Analysis Complete ===\n")
sink()
