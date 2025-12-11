# HCCbone: Spatial Transcriptomics and IMC Analysis of HCC Bone Metastasis

[![R](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://cran.r-project.org/)
[![SPATA2](https://img.shields.io/badge/SPATA2-Spatial%20Analysis-green.svg)](https://themilolab.github.io/SPATA2/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository contains a comprehensive spatial multi-omics analysis pipeline for characterizing the tumor microenvironment (TME) of **Hepatocellular Carcinoma (HCC) bone metastasis**. The analysis integrates:

- **Visium HD Spatial Transcriptomics** (8µm and 16µm resolution)
- **Imaging Mass Cytometry (IMC)** following the [Bodenmiller pipeline](https://bodenmillergroup.github.io/IMCDataAnalysis/)

The project was developed to support a case report investigating the unexpected ICI response in an HCC patient harboring a **CTNNB1 p.S23_S33del mutation**, which typically confers resistance to immunotherapy.

## Scientific Questions

This analysis addresses six key biological questions:

### 1. CD8⁺ T-cell Spatial Integration (IMC ↔ ST)
- Correlate IMC-derived CD8⁺ T-cell densities with CD8A/CD8B expression in spatial transcriptomics
- Validate CD8-rich regions across modalities

### 2. CD8⁺ T-cell Functional States
Characterize the functional landscape of tumor-infiltrating CD8⁺ T cells:

| State | Markers |
|-------|---------|
| **Effector/Cytotoxic** | IFNG, GZMB, PRF1 |
| **Exhaustion** | PDCD1 (PD-1), LAG3, HAVCR2 (TIM-3), TIGIT, TOX |
| **Tissue-Resident Memory (TRM)** | ITGAE (CD103), CD69, CXCR6 |

Key co-expression patterns: PD-1/CD69, PD-1/CD103

### 3. Myeloid and Stromal Niches
Analyze the spatial context surrounding CD8⁺ T-cell clusters:

| Population | Markers |
|------------|---------|
| **M1 Macrophages** | IL1B, NOS2, TNF, CD80, CD86 |
| **M2 Macrophages** | MRC1 (CD206), MSR1 (CD204), CD163 |
| **CAFs** | ACTA2 (α-SMA), POSTN, FAP, COL1A1 |

### 4. Salié Immunotype Features
Based on [Salié et al. (2023)](https://doi.org/10.1016/j.cell.2023.01.014):
- COL4A1 spatial relationship with CD8⁺ T cells
- DC neighborhood analysis (XCR1, CLEC9A)
- Tertiary Lymphoid Structure (TLS) detection: CXCL13, CCL19, MS4A1 (CD20), CR2, FCER2, CR1
- PD-1/CD20 co-localization

### 5. CTNNB1 Mutation Consequences
Investigate transcriptomic effects of the p.S23_S33del mutation:

| Pathway | Markers |
|---------|---------|
| **IFNγ Response** | STAT1, IRF1, ISG15 |
| **Antigen Presentation (HLA-I)** | HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2 |
| **Wnt/β-catenin Activity** | AXIN2, LEF1, TCF7, NKD1, RNF43, MYC, CCND1 |

### 6. NK Cell Neighborhood Analysis
- Identify NK-rich regions (NKG7, GNLY, KLRD1, KLRF1)
- Assess proximity to macrophages, CD8⁺/CD4⁺ T cells, mast cells
- Detect recurrent multicellular niche patterns

## Project Structure

```
HCCbone/
├── analysis/
│   ├── logs/                    # Execution logs
│   ├── objects/                 # Processed R objects
│   │   ├── *.rds               # SPATA2 objects
│   │   └── imc/                # IMC processed data
│   ├── plots/                   # All visualizations
│   │   ├── celltype/           # Cell type annotations
│   │   ├── cd8_spatial/        # CD8 spatial enrichment
│   │   ├── cd8_functional/     # CD8 functional states
│   │   ├── myeloid_caf/        # Myeloid & CAF niches
│   │   ├── salie/              # Salié immunotype features
│   │   ├── wnt_ifng/           # Wnt/IFNγ pathway analysis
│   │   ├── nk/                 # NK cell analysis
│   │   ├── st_spatial/         # ST spatial interactions
│   │   └── imc/                # IMC spatial analysis
│   └── tables/                  # CSV/TSV outputs
│       ├── global_summary/     # Cross-sample summaries
│       └── {analysis_type}/    # Per-analysis outputs
├── data/
│   └── VisiumHD/               # Raw Visium HD data
├── results/
│   └── spaceranger/            # SpaceRanger outputs
├── scripts/
│   ├── spata2/                 # Main analysis scripts
│   │   ├── plot_themes.R       # Color palettes & themes
│   │   ├── 00_*.R              # Data loading & QC
│   │   ├── 01-07_*.R           # ST analysis pipeline
│   │   ├── 07a_*.R             # ST spatial interactions
│   │   ├── 08-10_*.R           # IMC analysis pipeline
│   │   └── 11_*.R              # Summary report
│   ├── spaceranger/            # SpaceRanger scripts
│   └── archive/                # Legacy scripts
├── run_full_analysis.sh        # Main pipeline script
├── run_spatial_update.sh       # Spatial-only update script
└── .github/
    └── copilot-instructions.md # AI assistant instructions
```

## Analysis Pipeline

### Spatial Transcriptomics (SPATA2)

| Script | Description |
|--------|-------------|
| `00_load_build_objects.R` | Load SpaceRanger → SPATA2 objects |
| `00a_qc_normalization.R` | QC, filtering, log-normalization |
| `00b_dimred_clustering.R` | PCA, UMAP, t-SNE, K-means |
| `01_celltype_annotation.R` | 21 HCC-specific cell type annotations |
| `02_cd8_spatial_enrichment.R` | CD8 score spatial mapping |
| `03_cd8_functional_states.R` | Effector, exhaustion, TRM signatures |
| `04_myeloid_caf_niches.R` | M1/M2 polarization, CAF analysis |
| `05_salie_immunotype_markers.R` | Salié features, TLS detection |
| `06_ctnnb1_wnt_ifng_analysis.R` | Wnt vs IFNγ vs HLA-I analysis |
| `07_nk_neighbourhood_analysis.R` | NK cell spatial context |
| `07a_st_spatial_interactions.R` | Interaction graphs, communities, patches |

### IMC Analysis (Bodenmiller-style)

| Script | Description |
|--------|-------------|
| `08_imc_preprocessing.R` | TIFF loading, arcsinh normalization |
| `09_imc_phenotyping_clustering.R` | K-means phenotyping |
| `10_imc_spatial_analysis.R` | Comprehensive spatial analysis |
| `11_summary_report.R` | Multi-modal summary |

### Spatial Analysis Components

The pipeline includes comprehensive spatial analysis for both ST and IMC:

1. **Spatial Interaction Graphs** - Cell-cell connectivity networks
2. **Cellular Neighborhood Analysis** - Neighbor composition profiles
3. **Spatial Community Detection** - DBSCAN-based tissue regions
4. **Patch Detection** - Contiguous cell type clusters
5. **Interaction Enrichment** - Permutation-based enrichment/avoidance
6. **Spatial Context Analysis** - Density mapping, hotspot detection

## Installation

### Prerequisites

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", 
                   "igraph", "dbscan", "reshape2", "tiff"))

# Bioconductor
BiocManager::install(c("SPATA2", "ComplexHeatmap"))
```

### Clone Repository

```bash
git clone https://github.com/yourusername/HCCbone.git
cd HCCbone
```

## Usage

### Full Pipeline (HPC/Slurm)

```bash
# Submit full analysis
sbatch run_full_analysis.sh

# Monitor progress
tail -f analysis/logs/full_analysis_*.out
squeue -u $USER
```

### Individual Scripts

```bash
# Load R module (if on HPC)
module load R/4.3.0

# Run specific analysis
Rscript scripts/spata2/01_celltype_annotation.R
```

### Quick Start (Interactive)

```r
# Load a processed object
library(SPATA2)
obj <- readRDS("analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm_kmeans.rds")

# Check cell type distribution
table(obj@meta_obs$celltype)

# Generate spatial plot
source("scripts/spata2/plot_themes.R")
coords <- getCoordsDf(obj)
plot_df <- merge(coords, obj@meta_obs, by = "barcodes")
ggplot(plot_df, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = DISCRETE_COLORS) +
  coord_fixed() + HQ_THEME()
```

## Data Sources

### Visium HD
- Location: `Q:\MED2-A-SpatialOmics\CytAssist AGBB\projects\VisiumHD\VisiumHD 27038840`
- Resolution: 8µm (primary), 16µm (validation)
- Samples:
  - `D1_27038840_HCCboneMet_seg` - HCC bone metastasis
  - `SRout_27038840_HCC_tumour` - Primary HCC tumor

### IMC
- Location: `/data/projects/HCCbone/data/IMC/27038840_Hiemens/`
- 9 samples (4 Biopsy_Spine, 5 Resection)
- Panel: Immune and stromal markers

## Cell Type Annotation

21 cell types annotated using curated marker panels:

| Category | Cell Types |
|----------|------------|
| **Tumor/Epithelial** | Hepatocyte, HCC_tumour, Cholangiocyte |
| **T cells** | CD8_T, CD4_T, Treg, T_cell |
| **Myeloid** | Macrophage, M1_macro, M2_macro, Monocyte, DC, pDC |
| **Other Immune** | B_cell, Plasma, NK |
| **Stromal** | Fibroblast, CAF, Endothelial, LSEC, Stellate |

## Output Examples

### Publication-Quality Plots
All plots are generated at 300 DPI with:
- High-contrast color palettes (viridis, inferno, RdBu diverging)
- White backgrounds for print compatibility
- Consistent theming via `HQ_THEME()`

### Key Outputs

| Analysis | Plot | Table |
|----------|------|-------|
| Cell Types | `celltype_spatial_map.png` | `celltype_counts.csv` |
| CD8 Enrichment | `cd8_high_regions.png` | `cd8_stats.csv` |
| M1/M2 Ratio | `m1_m2_ratio.png` | `myeloid_caf_stats.csv` |
| TLS Detection | `tls_hotspots.png` | `tls_markers.csv` |
| Wnt/IFNγ | `wnt_ifng_comparison.png` | `pathway_correlations.csv` |
| Interactions | `interaction_enrichment.png` | `significant_interactions.csv` |

## HPC Configuration

Tested on Slurm cluster with these resource allocations:

| Partition | Time | CPUs | Memory | Use Case |
|-----------|------|------|--------|----------|
| short | 2h | ≤10 | 64G | Testing |
| long | 8h | 2-10 | 64-128G | Full analysis |

## References

1. **SPATA2**: [Kueckelhaus et al., 2023](https://themilolab.github.io/SPATA2/)
2. **Bodenmiller IMC Pipeline**: [bodenmillergroup.github.io/IMCDataAnalysis](https://bodenmillergroup.github.io/IMCDataAnalysis/)
3. **Salié Immunotypes**: [Salié et al., Cell 2023](https://doi.org/10.1016/j.cell.2023.01.014)
4. **CTNNB1 in HCC**: Review of β-catenin signaling in hepatocellular carcinoma

## Authors

- **Iman Sadeghi** - Spatial transcriptomics analysis
- **Maryam** - Case report presentation (GASL)
- **Ira** - IMC analysis, tumor vs metastasis comparison
- **Julia** - Project coordination, biological interpretation

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- AG Bengsch Lab for computational resources
- MED2-A-SpatialOmics for Visium HD data
- The Bodenmiller Lab for the IMC analysis framework

---

For questions or collaboration inquiries, please open an issue or contact the authors.
