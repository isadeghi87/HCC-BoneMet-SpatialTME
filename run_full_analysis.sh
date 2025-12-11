#!/bin/bash
#SBATCH --job-name=HCC_full_analysis
#SBATCH --output=analysis/logs/full_analysis_%j.out
#SBATCH --error=analysis/logs/full_analysis_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --partition=long
#SBATCH --qos=long

# HCC Bone Metastasis - Full Analysis Pipeline
# Runs all SPATA2 and IMC analysis scripts in sequence

set -e

cd /home/sadeghii/projects/HCCbone
mkdir -p analysis/logs

echo "=== HCC Full Analysis Pipeline ==="
echo "Start time: $(date)"

# Load R
module load R/4.3.0 2>/dev/null || true

SCRIPTS_DIR="scripts/spata2"

run_script() {
    local script=$1
    echo ""
    echo "=========================================="
    echo "Running: $script"
    echo "Time: $(date)"
    echo "=========================================="
    
    Rscript "$SCRIPTS_DIR/$script"
    
    if [ $? -eq 0 ]; then
        echo "✓ $script completed successfully"
    else
        echo "✗ $script failed"
        # Continue with next script instead of exiting
    fi
}

# ST Analysis Pipeline
echo ""
echo "=== SPATA2 Analysis Pipeline ==="

# Step 0: Data Loading and Preprocessing (skip if objects exist)
if [ ! -f "analysis/objects/D1_27038840_HCCboneMet_seg_raw.rds" ]; then
    run_script "00_load_build_objects.R"
fi

if [ ! -f "analysis/objects/D1_27038840_HCCboneMet_seg_norm_lognorm.rds" ]; then
    run_script "00a_qc_normalization.R"
fi

# Step 0b: Dimensionality Reduction (run if needed)
run_script "00b_dimred_clustering.R"

# Step 1: Cell Type Annotation
run_script "01_celltype_annotation.R"

# 1. CD8 Spatial Enrichment
run_script "02_cd8_spatial_enrichment.R"

# 2. CD8 Functional States
run_script "03_cd8_functional_states.R"

# 3. Myeloid and CAF Niches
run_script "04_myeloid_caf_niches.R"

# 4. Salié Immunotype Markers
run_script "05_salie_immunotype_markers.R"

# 5. CTNNB1/Wnt and IFNγ Analysis
run_script "06_ctnnb1_wnt_ifng_analysis.R"

# 6. NK Cell Neighbourhood Analysis
run_script "07_nk_neighbourhood_analysis.R"

# 7. Comprehensive ST Spatial Interactions (NEW)
run_script "07a_st_spatial_interactions.R"

# IMC Analysis Pipeline
echo ""
echo "=== IMC Analysis Pipeline ==="

# 7. IMC Preprocessing
run_script "08_imc_preprocessing.R"

# 8. IMC Phenotyping and Clustering
run_script "09_imc_phenotyping_clustering.R"

# 9. IMC Spatial Context Analysis
run_script "10_imc_spatial_analysis.R"

# Summary Report
echo ""
echo "=== Generating Summary Report ==="
run_script "11_summary_report.R"

echo ""
echo "=== Pipeline Complete ==="
echo "End time: $(date)"
echo ""
echo "Results:"
echo "  Tables: analysis/tables/"
echo "  Plots: analysis/plots/"
echo "  Summary: analysis/tables/global_summary/"
echo "  Logs: analysis/logs/"
