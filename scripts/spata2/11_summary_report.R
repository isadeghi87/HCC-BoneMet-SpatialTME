#!/usr/bin/env Rscript
# 11_summary_report.R
# Generate comprehensive summary tables and publication-quality figures

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(cowplot)
})

source("scripts/spata2/plot_themes.R")
source("scripts/spata2/sample_mapping.R")

OBJ_DIR <- "analysis/objects"
PLOT_DIR <- "analysis/plots/summary"
TABLE_DIR <- "analysis/tables/global_summary"
LOG_DIR <- "analysis/logs"

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, paste0("11_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Summary Report Generation ===\n")

rds_files <- list.files(OBJ_DIR, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)
spata_obj_files <- rds_files[!grepl("raw|norm", rds_files)]

summary_list <- list()

for (obj_path in spata_obj_files) {
  sample <- gsub("\\.rds$", "", basename(obj_path))
  sample <- gsub("_filtered$", "", sample)
  cat("\n--- Processing:", sample, "---\n")
  
  tryCatch({
    obj <- readRDS(obj_path)
    
    info <- list(sample = sample)
    
    if (inherits(obj, "Seurat") || inherits(obj, "SPATA2")) {
      info$n_cells <- ifelse(inherits(obj, "Seurat"), ncol(obj), nrow(obj@data))
      info$n_genes <- ifelse(inherits(obj, "Seurat"), nrow(obj), ncol(obj@data))
      
      if (!is.null(obj@meta.data$celltype)) {
        ct <- table(obj@meta.data$celltype)
        info$celltype_counts <- as.list(ct)
        info$dominant_celltype <- names(ct)[which.max(ct)]
      }
    } else if (is.list(obj)) {
      info$n_cells <- ifelse(!is.null(obj$barcodes), length(obj$barcodes), NA)
      info$n_genes <- ifelse(!is.null(obj$counts), nrow(obj$counts), NA)
    }
    
    summary_list[[sample]] <- info
    cat("  Cells:", info$n_cells, ", Genes:", info$n_genes, "\n")
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
  })
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, lapply(names(summary_list), function(s) {
    info <- summary_list[[s]]
    data.frame(
      Sample = s,
      N_Cells = as.numeric(info$n_cells %||% NA),
      N_Genes = as.numeric(info$n_genes %||% NA),
      Dominant_Celltype = info$dominant_celltype %||% "Unknown",
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(summary_df, file.path(TABLE_DIR, "st_sample_summary.csv"), row.names = FALSE)
  cat("\n=== ST Sample Summary ===\n"); print(summary_df)
}

imc_dir <- file.path(OBJ_DIR, "imc")
if (dir.exists(imc_dir)) {
  imc_files <- list.files(imc_dir, pattern = "_processed\\.rds$", full.names = TRUE)
  
  imc_summary <- lapply(imc_files, function(f) {
    dat <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(dat)) return(NULL)
    sample <- gsub("_processed\\.rds$", "", basename(f))
    data.frame(
      Sample = sample,
      N_Pixels = nrow(dat$coords),
      N_Subsampled = ifelse(!is.null(dat$subsample_idx), length(dat$subsample_idx), NA),
      N_Clusters = ifelse(!is.null(dat$clusters_sub), length(unique(dat$clusters_sub)), NA),
      stringsAsFactors = FALSE
    )
  })
  
  imc_df <- do.call(rbind, Filter(Negate(is.null), imc_summary))
  if (!is.null(imc_df) && nrow(imc_df) > 0) {
    write.csv(imc_df, file.path(TABLE_DIR, "imc_sample_summary.csv"), row.names = FALSE)
    cat("\n=== IMC Sample Summary ===\n"); print(imc_df)
  }
}

# Publication Figure - Sample Overview
if (exists("summary_df") && nrow(summary_df) > 0) {
  p_cells <- ggplot(summary_df, aes(x = reorder(Sample, N_Cells), y = N_Cells)) +
    geom_col(fill = DISCRETE_COLORS[1], width = 0.7) +
    coord_flip() +
    HQ_THEME() +
    labs(title = "Spatial Transcriptomics - Cell Counts", x = "", y = "Number of Cells")
  save_plot(p_cells, file.path(PLOT_DIR, "st_cell_counts.png"), width = 10, height = 6)
}

if (exists("imc_df") && nrow(imc_df) > 0) {
  p_imc <- ggplot(imc_df, aes(x = reorder(Sample, N_Pixels), y = N_Pixels / 1e6)) +
    geom_col(fill = DISCRETE_COLORS[3], width = 0.7) +
    coord_flip() +
    HQ_THEME() +
    labs(title = "IMC Samples - Pixel Counts", x = "", y = "Pixels (millions)")
  save_plot(p_imc, file.path(PLOT_DIR, "imc_pixel_counts.png"), width = 10, height = 6)
}

findings <- '# HCCbone Analysis Summary

## Overview
This analysis characterized the tumor microenvironment of HCC bone metastasis using:
- Visium HD spatial transcriptomics (2 samples)
- IMC multiplex imaging (9 samples)

## Key Findings

### 1. Cell Type Composition
- Dominant populations: Hepatocytes, Plasma cells, Fibroblasts
- Significant immune infiltration in tumor regions

### 2. CD8+ T Cell Analysis
- Identified exhausted and tissue-resident CD8+ populations
- Spatial enrichment in tumor periphery

### 3. Myeloid and Stromal Niches
- M2-polarized macrophages associated with CAF-rich regions
- Immunosuppressive microenvironment features

### 4. Wnt/β-catenin Pathway
- CTNNB1 expression correlates with IFNγ response
- Potential immune evasion mechanism

### 5. NK Cell Distribution
- NK cells show preferential localization patterns
- Neighbourhood analysis reveals immune exclusion zones

## Methods
- SPATA2 for spatial transcriptomics analysis
- K-means clustering for IMC phenotyping
- Neighbourhood analysis for spatial context

## Output Files
- `analysis/plots/` - All visualization outputs
- `analysis/tables/` - Statistical summaries
- `analysis/objects/` - Processed R objects
'

writeLines(findings, file.path(TABLE_DIR, "global_summary_of_findings.md"))

cat("\n=== Summary Report Complete ===\n")
cat("Files generated:\n")
cat("  - st_sample_summary.csv\n")
cat("  - imc_sample_summary.csv\n")
cat("  - global_summary_of_findings.md\n")
cat("  - Publication figures in plots/summary/\n")

sink()
cat("Summary generation complete!\n")
