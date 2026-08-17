#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=Master_QC
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/Master_QC_%j.out
#SBATCH --error=/path/to/project_directory/logs/Master_QC_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment for MultiQC
source /path/to/virtualenv/bin/activate

FASTQC="/path/to/software/FastQC/fastqc"
BASE_DIR="/path/to/project_directory"
QC_DIR="${BASE_DIR}/QC_Reports"

# Ensure output directory exists
mkdir -p ${QC_DIR}

# --- 2. EXECUTE FASTQC ---
echo "========================================="
echo "STEP 1: Running FastQC on all Merged Files"
echo "========================================="

# Find all merged FASTQ files across patient folders and process concurrently
find ${BASE_DIR}/*_analysis -name "*_merged.fastq.gz" | xargs -I {} -P 8 ${FASTQC} {} -o ${QC_DIR}

# --- 3. AGGREGATE WITH MULTIQC ---
echo "========================================="
echo "STEP 2: Generating MultiQC Report"
echo "========================================="

cd ${QC_DIR}
multiqc .

echo "Quality Control Complete. MultiQC report generated successfully."
