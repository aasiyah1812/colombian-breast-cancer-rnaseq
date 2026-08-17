#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=GenomicsDB_Joint_Genotyping
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/GenomicsDB_58_%j.out
#SBATCH --error=/path/to/project_directory/logs/GenomicsDB_58_%j.err
#SBATCH --time=500:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment
source /path/to/virtualenv/bin/activate

BASE_DIR="/path/to/project_directory"
GVCF_DIR="${BASE_DIR}/All_GVCFs"
REF="${BASE_DIR}/GRCh38.p14.genome.fa"

OUTPUT_DIR="${BASE_DIR}/Joint_Variant_Calling_GenomicsDB"
mkdir -p ${OUTPUT_DIR}

# Create a dedicated temporary directory to manage heavy I/O memory load
TMP_DIR="${BASE_DIR}/tmp_gatk"
mkdir -p ${TMP_DIR}

DB_WORKSPACE="${OUTPUT_DIR}/full_58_workspace"
FINAL_VCF="${OUTPUT_DIR}/full_58_GenomicsDB_variants.vcf.gz"

# Clear any previous workspace to prevent database corruption
rm -rf ${DB_WORKSPACE}

echo "========================================="
echo "STARTING JOINT GENOTYPING FOR 58 PATIENTS"
echo "========================================="

# --- 2. BUILD INPUT LIST AUTOMATICALLY ---
# Dynamically collect all individual GVCFs into a single input string
INPUT_ARGS=""
for file in ${GVCF_DIR}/*.g.vcf.gz; do
    INPUT_ARGS="$INPUT_ARGS -V $file"
done

# --- 3. BUILD CONTIG FILTER STRING ---
# Restrict analysis to main chromosomes to optimize computational efficiency
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
MAIN_CONTIGS=""
for chrom in $CHROMS; do
  MAIN_CONTIGS="$MAIN_CONTIGS -L $chrom"
done

# --- 4. RUN GenomicsDBImport ---
echo "Step 1: Importing GVCFs to GenomicsDB..."
# Utilizing explicit Java memory allocation and batching to manage high memory demands
gatk --java-options "-Xmx60g" GenomicsDBImport \
    --genomicsdb-workspace-path ${DB_WORKSPACE} \
    --tmp-dir ${TMP_DIR} \
    ${INPUT_ARGS} \
    ${MAIN_CONTIGS} \
    --batch-size 25 
    
# --- 5. GENOTYPE GVCFs ---
echo "Step 2: Performing Joint Genotyping..."
gatk --java-options "-Xmx60g" GenotypeGVCFs \
    -R ${REF} \
    -V gendb://${DB_WORKSPACE} \
    -O ${FINAL_VCF} \
    --tmp-dir ${TMP_DIR} \
    --standard-min-confidence-threshold-for-calling 20

# Clean up the temporary directory
rm -rf ${TMP_DIR}

echo "Joint Genotyping Complete. Final cohort VCF generated successfully."
