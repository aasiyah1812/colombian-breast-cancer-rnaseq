#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=Merge_FASTQ_014B
#SBATCH --partition=fuchs
#SBATCH --output=014B_analysis/merge_014B.out
#SBATCH --error=014B_analysis/merge_014B.err
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

echo "Starting FASTQ merge for Patient 014B"
date

# --- Script Content ---
PATIENT="014B" 
RAW_DATA_DIR="./RAW_FASTQ_DATA"
OUTPUT_DIR="${PATIENT}_analysis"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# 1. Merge Read 1 (R1) files across all 4 sequencing lanes
echo "Merging R1 files (L001-L004)..."
zcat ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L001_R1_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L002_R1_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L003_R1_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L004_R1_001.fastq.gz | gzip > ${OUTPUT_DIR}/${PATIENT}_R1_merged.fastq.gz

# 2. Merge Read 2 (R2) files across all 4 sequencing lanes
echo "Merging R2 files (L001-L004)..."
zcat ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L001_R2_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L002_R2_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L003_R2_001.fastq.gz \
     ${RAW_DATA_DIR}/${PATIENT}-Tn_S14_L004_R2_001.fastq.gz | gzip > ${OUTPUT_DIR}/${PATIENT}_R2_merged.fastq.gz

echo "Merge complete for ${PATIENT}. Checking file sizes."
ls -lh ${OUTPUT_DIR}/${PATIENT}_R*_merged.fastq.gz
date
