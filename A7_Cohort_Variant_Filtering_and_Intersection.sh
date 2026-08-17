#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=Filter_Intersect_VCF
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/Filter_Intersect_%j.out
#SBATCH --error=/path/to/project_directory/logs/Filter_Intersect_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment
source /path/to/virtualenv/bin/activate

BASE_DIR="/path/to/project_directory"
COHORT_VCF="${BASE_DIR}/01_Cohort_Data/full_58_GenomicsDB_variants.vcf.gz"
REF_1000G="${BASE_DIR}/02_Reference_Data/True_1000G_Reference.vcf.gz"
OUT_DIR="${BASE_DIR}/05_Pipeline_Results"

BCFTOOLS="/path/to/software/bcftools"

# Ensure output directory exists
mkdir -p ${OUT_DIR}

# --- 2. HARD FILTERING ---
echo "Step 1: Applying quality filters and isolating biallelic SNPs..."
# Apply standard hard filters for quality and restrict to standard biallelic SNPs
${BCFTOOLS} filter -i 'QD>2.0 && FS<60.0 && MQ>40.0 && SOR<3.0' ${COHORT_VCF} | \
${BCFTOOLS} view -m2 -M2 -v snps -O z -o ${OUT_DIR}/cohort_58_filtered.vcf.gz

# Index the filtered file for downstream processing
${BCFTOOLS} index ${OUT_DIR}/cohort_58_filtered.vcf.gz

# --- 3. GERMLINE INTERSECTION ---
echo "Step 2: Intersecting cohort variants with the 1000 Genomes reference..."
# Retain exclusively those variants present in the 1000G reference database
# This isolates true germline ancestry markers from somatic mutations
${BCFTOOLS} isec -n =2 -w 1 \
    ${OUT_DIR}/cohort_58_filtered.vcf.gz \
    ${REF_1000G} \
    -O z -o ${OUT_DIR}/cohort_58_intersected.vcf.gz

echo "Variant filtering and intersection complete. Outputs saved successfully."
