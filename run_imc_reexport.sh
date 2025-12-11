#!/bin/bash
#SBATCH --job-name=imc_reexp
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=analysis/logs/imc_reexport_%j.out
#SBATCH --error=analysis/logs/imc_reexport_%j.err

# Re-export IMC outputs with standardized naming

cd /home/sadeghii/projects/HCCbone

echo "=== IMC Re-export with Standardized Naming ==="
echo "Time: $(date)"
echo ""

module load R/4.3.0 2>/dev/null || true

Rscript scripts/spata2/99b_reexport_imc_outputs.R

echo ""
echo "=== IMC Re-export Complete ==="
echo "Time: $(date)"
