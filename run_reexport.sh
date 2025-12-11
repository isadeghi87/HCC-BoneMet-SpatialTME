#!/bin/bash
#SBATCH --job-name=reexport
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=analysis/logs/reexport_%j.out
#SBATCH --error=analysis/logs/reexport_%j.err

# Re-export outputs with standardized naming
# Fast job - just loads objects and regenerates plots/tables

cd /home/sadeghii/projects/HCCbone

echo "=== Re-export with Standardized Naming ==="
echo "Time: $(date)"
echo "Node: $(hostname)"
echo ""

# Load R
module load R/4.3.0 2>/dev/null || true

# Run ST re-export script
echo "=== Part 1: ST Re-export ==="
Rscript scripts/spata2/99_reexport_standardized_outputs.R

# Run IMC re-export script
echo ""
echo "=== Part 2: IMC Re-export ==="
Rscript scripts/spata2/99b_reexport_imc_outputs.R

echo ""
echo "=== All Re-exports Complete ==="
echo "Time: $(date)"
