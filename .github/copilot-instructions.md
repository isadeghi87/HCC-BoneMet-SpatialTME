# Copilot Instructions for HCCbone Project

## Project Overview

This is a **Hepatocellular Carcinoma (HCC) Bone Metastasis** spatial transcriptomics and imaging mass cytometry analysis project. The project uses SPATA2 for Visium HD analysis and custom R pipelines for IMC data.

## Data Structure

- **Visium HD data**: `/home/sadeghii/projects/HCCbone/results/spaceranger/`
- **IMC data**: `/data/projects/HCCbone/data/IMC/27038840_Hiemens/`
- **Analysis outputs**: `analysis/plots/`, `analysis/tables/`, `analysis/objects/`

## Key Analysis Components

### Spatial Transcriptomics (ST)
1. **Cell Type Annotation**: 21 HCC-specific cell types using marker panels
2. **CD8+ T Cell Analysis**: Spatial enrichment, exhaustion markers
3. **Myeloid/Stromal Analysis**: TAMs, M1/M2 polarization, CAFs
4. **Pathway Analysis**: Wnt/β-catenin, IFNγ response
5. **NK Cell Analysis**: Neighbourhood analysis

### Imaging Mass Cytometry (IMC)
1. **Preprocessing**: TIFF loading, arcsinh normalization
2. **Phenotyping**: K-means clustering
3. **Spatial Analysis**: Neighbourhood matrices, interaction testing

## Coding Standards

### R Scripts
- Use `suppressPackageStartupMessages()` for clean loading
- Source `scripts/spata2/plot_themes.R` for consistent plotting
- Use `HQ_THEME()` and high-contrast color palettes
- Save plots at **DPI 300** minimum for publication quality
- Log all outputs to `analysis/logs/`

### Plot Quality
```r
# Use these color palettes
DISCRETE_COLORS  # 20 high-contrast colors
GRADIENT_INFERNO # Expression heatmaps
GRADIENT_VIRIDIS # General gradients
DIVERGENT_RDBU   # Correlation/change plots

# Standard theme
HQ_THEME()  # Publication-ready base theme

# Save function
save_plot(plot, filepath, width = 12, height = 10, dpi = 300)
```

### HPC Execution
- **Short partition**: 2 hours, 10 CPUs max - for testing
- **Medium partition**: 48 hours, 20 CPUs max - for full analysis
- **Long partition**: 7 days, 30 CPUs max - for large jobs

## File Naming Conventions
- Analysis scripts: `NN_descriptive_name.R` (e.g., `01_celltype_annotation.R`)
- Plots: `sample_analysis_type.png`
- Tables: `sample_data_type.csv`

## Key Cell Type Markers

### Immune Cells
- CD8 T cells: CD8A, CD8B, CD3D, GZMB
- NK cells: NCAM1, NKG7, GNLY
- Macrophages: CD68, CD163, CSF1R
- Tregs: FOXP3, IL2RA

### Stromal Cells
- CAFs: ACTA2, FAP, PDGFRA
- Endothelial: PECAM1, CDH5

### Tumor Cells
- Hepatocytes: ALB, APOA1, CYP3A4
- Cholangiocytes: KRT19, EPCAM

## Common Tasks

### Load SPATA2 Object
```r
library(SPATA2)
obj <- readRDS("analysis/objects/D1_27038840_HCCboneMet_seg_filtered.rds")
```

### Generate Expression Heatmap
```r
source("scripts/spata2/plot_themes.R")
p <- ggplot(data, aes(x, y, fill = expression)) +
  geom_tile() +
  scale_fill_gradientn(colors = GRADIENT_INFERNO) +
  HQ_THEME()
save_plot(p, "output.png")
```

### Submit HPC Job
```bash
sbatch run_full_analysis.sh
squeue -u $USER  # Check status
scancel <job_id>  # Cancel if needed
```

## Dependencies
- R 4.3.0+
- SPATA2
- ggplot2, dplyr, tidyr
- Seurat (for some analyses)
- tiff (for IMC)

## Contact
Analysis pipeline maintained for HCC bone metastasis TME characterization study.
