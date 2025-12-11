#!/usr/bin/env Rscript
# 07a_st_spatial_interactions.R
# Comprehensive ST Spatial Analysis: Interaction graphs, communities, patches
# HIGH-QUALITY PLOTTING VERSION

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(igraph)
  library(dbscan)
  library(patchwork)
})

source("scripts/spata2/plot_themes.R")

PLOT_DIR <- "analysis/plots/st_spatial"
TABLE_DIR <- "analysis/tables/st_spatial"
LOG_DIR <- "analysis/logs"

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("07a_st_spatial_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Comprehensive ST Spatial Analysis ===\n")
cat("Components: Interaction graphs, Spatial communities, Patch detection,\n")
cat("            Cellular neighborhoods, Interaction enrichment\n\n")

# ============================================================================
# 1. SPATIAL INTERACTION GRAPH
# ============================================================================

build_st_spatial_graph <- function(coords, radius = 50) {
  cat("  Building spatial interaction graph (radius =", radius, "µm)...\n")
  n <- nrow(coords)
  edges <- list()
  
  sample_size <- min(n, 20000)
  sample_idx <- sample(n, sample_size)
  
  for (i in sample_idx) {
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    for (j in neighbors) {
      edges[[length(edges) + 1]] <- c(i, j)
    }
  }
  
  if (length(edges) == 0) return(NULL)
  
  edge_matrix <- do.call(rbind, edges)
  g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)
  g <- igraph::simplify(g)
  
  cat("    Nodes:", igraph::vcount(g), "Edges:", igraph::ecount(g), "\n")
  return(g)
}

# ============================================================================
# 2. CELLULAR NEIGHBORHOOD COMPOSITION
# ============================================================================

compute_st_neighborhoods <- function(coords, celltypes, radius = 50) {
  cat("  Computing cellular neighborhoods...\n")
  n <- nrow(coords)
  unique_types <- sort(unique(celltypes))
  n_types <- length(unique_types)
  
  cn_matrix <- matrix(0, nrow = n, ncol = n_types, dimnames = list(NULL, unique_types))
  
  sample_idx <- if (n > 20000) sample(n, 20000) else seq_len(n)
  
  for (idx in seq_along(sample_idx)) {
    i <- sample_idx[idx]
    if (idx %% 5000 == 0) cat("    Processed", idx, "/", length(sample_idx), "\n")
    
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    
    if (length(neighbors) > 0) {
      neighbor_types <- celltypes[neighbors]
      for (ct in unique_types) {
        cn_matrix[i, ct] <- sum(neighbor_types == ct)
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

detect_st_communities <- function(coords, celltypes, eps = 80, minPts = 10) {
  cat("  Detecting spatial communities (DBSCAN)...\n")
  
  db_result <- dbscan::dbscan(as.matrix(coords[, c("x", "y")]), eps = eps, minPts = minPts)
  communities <- db_result$cluster
  n_communities <- max(communities)
  
  cat("    Found", n_communities, "spatial communities\n")
  cat("    Noise points:", sum(communities == 0), "\n")
  
  # Community composition
  community_composition <- data.frame()
  for (comm in 1:n_communities) {
    comm_cells <- which(communities == comm)
    if (length(comm_cells) > 0) {
      comm_types <- celltypes[comm_cells]
      comp <- table(comm_types) / length(comm_cells)
      comp_df <- data.frame(community = comm, celltype = names(comp), proportion = as.numeric(comp))
      community_composition <- rbind(community_composition, comp_df)
    }
  }
  
  return(list(membership = communities, n_communities = n_communities, composition = community_composition))
}

# ============================================================================
# 4. PATCH DETECTION
# ============================================================================

detect_st_patches <- function(coords, celltypes, eps = 40, minPts = 5) {
  cat("  Detecting cell type patches...\n")
  
  unique_types <- sort(unique(celltypes))
  patch_results <- list()
  
  for (ct in unique_types) {
    ct_idx <- which(celltypes == ct)
    if (length(ct_idx) < minPts) next
    
    ct_coords <- coords[ct_idx, c("x", "y")]
    
    db <- tryCatch({
      dbscan::dbscan(as.matrix(ct_coords), eps = eps, minPts = minPts)
    }, error = function(e) NULL)
    
    if (!is.null(db)) {
      n_patches <- max(db$cluster)
      patch_sizes <- table(db$cluster[db$cluster > 0])
      
      patch_results[[ct]] <- list(
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
# 5. INTERACTION ENRICHMENT (Permutation test)
# ============================================================================

compute_st_interaction_enrichment <- function(coords, celltypes, radius = 50, n_permutations = 50) {
  cat("  Computing interaction enrichment...\n")
  
  unique_types <- sort(unique(celltypes))
  n_types <- length(unique_types)
  n <- nrow(coords)
  
  # Observed
  observed <- matrix(0, nrow = n_types, ncol = n_types, dimnames = list(unique_types, unique_types))
  
  sample_idx <- if (n > 15000) sample(n, 15000) else seq_len(n)
  
  for (i in sample_idx) {
    dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
    neighbors <- which(dists > 0 & dists <= radius)
    ci <- celltypes[i]
    for (j in neighbors) {
      cj <- celltypes[j]
      observed[ci, cj] <- observed[ci, cj] + 1
    }
  }
  
  # Permutations
  cat("    Running", n_permutations, "permutations...\n")
  perm_results <- array(0, dim = c(n_types, n_types, n_permutations))
  
  for (p in 1:n_permutations) {
    if (p %% 10 == 0) cat("      Permutation", p, "\n")
    
    perm_types <- sample(celltypes)
    perm_mat <- matrix(0, nrow = n_types, ncol = n_types)
    
    for (i in sample_idx[1:min(8000, length(sample_idx))]) {
      dists <- sqrt((coords$x - coords$x[i])^2 + (coords$y - coords$y[i])^2)
      neighbors <- which(dists > 0 & dists <= radius)
      ci <- perm_types[i]
      for (j in neighbors) {
        cj <- perm_types[j]
        ci_idx <- which(unique_types == ci)
        cj_idx <- which(unique_types == cj)
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
  
  dimnames(z_scores) <- list(unique_types, unique_types)
  dimnames(p_values) <- list(unique_types, unique_types)
  
  return(list(observed = observed, expected = perm_means, z_scores = z_scores, p_values = p_values))
}

# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

analyze_st_spatial <- function(obj_path, sample_name) {
  cat("\n", strrep("=", 60), "\n")
  cat("Processing:", sample_name, "\n")
  cat(strrep("=", 60), "\n")
  
  obj <- readRDS(obj_path)
  coords <- getCoordsDf(obj)
  meta <- obj@meta_obs
  
  if (!"celltype" %in% colnames(meta)) {
    cat("No celltype annotation, skipping spatial analysis.\n")
    return(NULL)
  }
  
  celltypes <- meta$celltype[match(coords$barcodes, meta$barcodes)]
  
  # Filter valid cells
  valid <- !is.na(celltypes) & celltypes != "Unassigned"
  coords <- coords[valid, ]
  celltypes <- celltypes[valid]
  
  cat("Cells for spatial analysis:", nrow(coords), "\n")
  cat("Cell types:", length(unique(celltypes)), "\n\n")
  
  sample_plot_dir <- file.path(PLOT_DIR, sample_name)
  sample_table_dir <- file.path(TABLE_DIR, sample_name)
  dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  # 1. Spatial Interaction Graph
  cat("\n[1/5] Spatial Interaction Graph\n")
  tryCatch({
    graph <- build_st_spatial_graph(coords, radius = 50)
    results$graph <- graph
    
    if (!is.null(graph)) {
      graph_stats <- data.frame(
        metric = c("nodes", "edges", "density", "avg_degree", "components"),
        value = c(igraph::vcount(graph), igraph::ecount(graph), igraph::edge_density(graph),
                  mean(igraph::degree(graph)), igraph::components(graph)$no)
      )
      write.csv(graph_stats, file.path(sample_table_dir, "graph_statistics.csv"), row.names = FALSE)
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 2. Cellular Neighborhoods
  cat("\n[2/5] Cellular Neighborhoods\n")
  tryCatch({
    cn <- compute_st_neighborhoods(coords, celltypes, radius = 50)
    results$neighborhoods <- cn
    
    cn_df <- as.data.frame(cn$normalized[cn$sample_idx, ])
    cn_df$celltype <- celltypes[cn$sample_idx]
    write.csv(cn_df, file.path(sample_table_dir, "cellular_neighborhoods.csv"), row.names = FALSE)
    
    # Summarize by cell type
    cn_summary <- cn_df %>%
      group_by(celltype) %>%
      summarise(across(everything(), mean, na.rm = TRUE)) %>%
      pivot_longer(-celltype, names_to = "neighbor_type", values_to = "proportion")
    
    # Filter to top 15 cell types
    top_types <- cn_df %>% count(celltype, sort = TRUE) %>% head(15) %>% pull(celltype)
    cn_summary_filt <- cn_summary %>% filter(celltype %in% top_types, neighbor_type %in% top_types)
    
    p_cn <- ggplot(cn_summary_filt, aes(x = celltype, y = neighbor_type, fill = proportion)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradientn(colors = GRADIENT_VIRIDIS) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = element_text(size = 9)) +
      labs(title = paste(sample_name, "- Cellular Neighborhood Composition"),
           x = "Cell Type", y = "Neighbor Type", fill = "Proportion")
    save_plot(p_cn, file.path(sample_plot_dir, "cellular_neighborhoods.png"), width = 12, height = 10)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 3. Spatial Communities
  cat("\n[3/5] Spatial Communities\n")
  tryCatch({
    communities <- detect_st_communities(coords, celltypes, eps = 80, minPts = 10)
    results$communities <- communities
    
    comm_df <- data.frame(x = coords$x, y = coords$y, celltype = celltypes, community = communities$membership)
    write.csv(comm_df, file.path(sample_table_dir, "spatial_communities.csv"), row.names = FALSE)
    
    if (nrow(communities$composition) > 0) {
      write.csv(communities$composition, file.path(sample_table_dir, "community_composition.csv"), row.names = FALSE)
    }
    
    # Plot communities
    plot_df <- comm_df %>% filter(community > 0)
    if (nrow(plot_df) > 50000) plot_df <- plot_df[sample(nrow(plot_df), 50000), ]
    
    if (nrow(plot_df) > 0) {
      p_comm <- ggplot(plot_df, aes(x = x, y = y, color = factor(community))) +
        geom_point(size = 0.2, alpha = 0.6) +
        scale_color_manual(values = DISCRETE_COLORS) +
        HQ_THEME() +
        coord_fixed() +
        labs(title = paste(sample_name, "- Spatial Communities"), x = "X (µm)", y = "Y (µm)", color = "Community") +
        guides(color = guide_legend(override.aes = list(size = 3)))
      save_plot(p_comm, file.path(sample_plot_dir, "spatial_communities.png"), width = 14, height = 12)
    }
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 4. Patch Detection
  cat("\n[4/5] Patch Detection\n")
  tryCatch({
    patches <- detect_st_patches(coords, celltypes, eps = 40, minPts = 5)
    results$patches <- patches
    
    patch_summary <- do.call(rbind, lapply(names(patches), function(ct) {
      data.frame(celltype = ct, n_patches = patches[[ct]]$n_patches,
                 mean_size = patches[[ct]]$mean_size, max_size = patches[[ct]]$max_size)
    }))
    write.csv(patch_summary, file.path(sample_table_dir, "patch_summary.csv"), row.names = FALSE)
    
    # Filter to top patches
    patch_summary_top <- patch_summary %>% filter(n_patches > 0) %>% arrange(desc(n_patches)) %>% head(20)
    
    p_patches <- ggplot(patch_summary_top, aes(x = reorder(celltype, -n_patches), y = n_patches, fill = mean_size)) +
      geom_col(width = 0.7) +
      scale_fill_gradientn(colors = GRADIENT_INFERNO) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(sample_name, "- Cell Type Patches"), x = "Cell Type", y = "Number of Patches", fill = "Mean Size")
    save_plot(p_patches, file.path(sample_plot_dir, "patch_detection.png"), width = 12, height = 7)
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  # 5. Interaction Enrichment
  cat("\n[5/5] Interaction Enrichment\n")
  tryCatch({
    # Use top cell types only for efficiency
    type_counts <- table(celltypes)
    top_types <- names(sort(type_counts, decreasing = TRUE))[1:min(12, length(type_counts))]
    
    top_idx <- celltypes %in% top_types
    coords_top <- coords[top_idx, ]
    celltypes_top <- celltypes[top_idx]
    
    interactions <- compute_st_interaction_enrichment(coords_top, celltypes_top, radius = 50, n_permutations = 30)
    results$interactions <- interactions
    
    write.csv(interactions$observed, file.path(sample_table_dir, "interaction_observed.csv"))
    write.csv(interactions$z_scores, file.path(sample_table_dir, "interaction_zscore.csv"))
    write.csv(interactions$p_values, file.path(sample_table_dir, "interaction_pvalues.csv"))
    
    # Z-score heatmap
    z_melt <- melt(interactions$z_scores)
    colnames(z_melt) <- c("Source", "Target", "Z_score")
    z_melt$Z_score <- pmax(pmin(z_melt$Z_score, 5), -5)
    
    p_zscore <- ggplot(z_melt, aes(x = Target, y = Source, fill = Z_score)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradientn(colors = DIVERGENT_RDBU, limits = c(-5, 5)) +
      HQ_THEME() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(sample_name, "- Cell Type Interaction Enrichment"),
           subtitle = "Z-score: >0 enriched (red), <0 avoided (blue)",
           x = "Target", y = "Source", fill = "Z-score")
    save_plot(p_zscore, file.path(sample_plot_dir, "interaction_enrichment.png"), width = 11, height = 9)
    
    # Significant interactions
    sig_int <- z_melt %>%
      left_join(melt(interactions$p_values) %>% rename(Source = Var1, Target = Var2, p_value = value), by = c("Source", "Target")) %>%
      filter(p_value < 0.05) %>%
      arrange(desc(abs(Z_score)))
    write.csv(sig_int, file.path(sample_table_dir, "significant_interactions.csv"), row.names = FALSE)
    cat("    Significant interactions:", nrow(sig_int), "\n")
  }, error = function(e) cat("  Error:", e$message, "\n"))
  
  cat("\n✓ ST spatial analysis complete for", sample_name, "\n")
  return(results)
}

# ============================================================================
# PROCESS SAMPLES
# ============================================================================

samples <- list(
  list(path = "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds", name = "D1_27038840_HCCboneMet_seg"),
  list(path = "analysis/objects/SRout_27038840_HCC_tumour_norm_lognorm_kmeans.rds", name = "SRout_27038840_HCC_tumour")
)

for (s in samples) {
  if (file.exists(s$path)) {
    tryCatch(analyze_st_spatial(s$path, s$name), error = function(e) cat("Error:", e$message, "\n"))
  } else {
    cat("Object not found:", s$path, "\n")
  }
}

cat("\n=== ST Spatial Analysis Complete ===\n")
sink()
