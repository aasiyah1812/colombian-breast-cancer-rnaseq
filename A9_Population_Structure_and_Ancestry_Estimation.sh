#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=ADMIXTURE_Array
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/admix_K%a_%j.out
#SBATCH --error=/path/to/project_directory/logs/admix_K%a_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=2-10

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment
source /path/to/virtualenv/bin/activate

BASE_DIR="/path/to/project_directory"
OUT_DIR="${BASE_DIR}/05_Pipeline_Results"
ADMIXTURE="/path/to/software/admixture_linux-1.3.0/admixture"

# Navigate to the output directory containing the final LD-pruned dataset
cd ${OUT_DIR}

echo "========================================="
echo "Running ADMIXTURE for K = $SLURM_ARRAY_TASK_ID"
echo "========================================="

# --- 2. EXECUTE ADMIXTURE ---
# Run ADMIXTURE utilizing 4 threads (-j4) and 10-fold cross-validation (--cv=10)
# Output is dynamically piped to a log file to extract CV error values downstream
${ADMIXTURE} --cv=10 FINAL_Ancestry_Data.bed $SLURM_ARRAY_TASK_ID -j4 | tee admixture_K${SLURM_ARRAY_TASK_ID}_log.txt

echo "ADMIXTURE modeling completed for K = $SLURM_ARRAY_TASK_ID"
