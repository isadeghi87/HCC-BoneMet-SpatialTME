#!/bin/bash
#SBATCH --job-name=HCC_spatial_update
#SBATCH --output=analysis/logs/spatial_update_%j.out
#SBATCH --error=analysis/logs/spatial_update_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --partition=long
#SBATCH --qos=long

# HCC Bone Metastasis - Spatial Analysis Update
# Runs only the remaining/updated spatial analysis scripts

set -e

cd /home/sadeghii/projects/HCCbone
mkdir -p analysis/logs

echo "=== HCC Spatial Analysis Update ==="
echo "Start time: $(date)"

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
    fi
}

# Skip already completed steps (01-06 completed in Job 133)
# Continue from NK analysis

echo ""
echo "=== Continuing from NK Analysis ==="

# 7. NK Neighbourhood Analysis
run_script "07_nk_neighbourhood_analysis.R"

# 7a. NEW: Comprehensive ST Spatial Interactions
run_script "07a_st_spatial_interactions.R"

# IMC Analysis Pipeline
echo ""
echo "=== IMC Analysis Pipeline ==="

# 8. IMC Preprocessing
run_script "08_imc_preprocessing.R"

# 9. IMC Phenotyping and Clustering
run_script "09_imc_phenotyping_clustering.R"

# 10. IMC Comprehensive Spatial Analysis (UPDATED)
run_script "10_imc_spatial_analysis.R"

# Summary Report
echo ""
echo "=== Generating Summary Report ==="
run_script "11_summary_report.R"

echo ""
echo "=== Spatial Update Complete ==="
echo "End time: $(date)"
echo ""
echo "Results:"
echo "  ST Spatial: analysis/plots/st_spatial/, analysis/tables/st_spatial/"
echo "  IMC Spatial: analysis/plots/imc/, analysis/tables/imc/"
echo "  Summary: analysis/tables/global_summary/"
