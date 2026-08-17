#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=GATK_HaplotypeCaller_GVCF
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/GVCF_%A_%a.out
#SBATCH --error=/path/to/project_directory/logs/GVCF_%A_%a.err
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --array=1-58%10

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment
source /path/to/virtualenv/bin/activate

BASE_DIR="/path/to/project_directory"

# Extract patient ID from the master list using the array task ID
PATIENT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BASE_DIR}/patient_list.txt)

if [ -z "$PATIENT" ]; then
    echo "Error: No patient found for Task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Set up specific patient output and reference paths
OUTPUT_DIR="${BASE_DIR}/${PATIENT}_analysis"
REF="${BASE_DIR}/GRCh38.p14.genome.fa"

# Centralized directory for all output GVCFs
GVCF_DIR="${BASE_DIR}/All_GVCFs"
mkdir -p $GVCF_DIR

# Define Input and Intermediate Files
STAR_BAM="${OUTPUT_DIR}/${PATIENT}_star_Aligned.sortedByCoord.out.bam"
RG_BAM="${OUTPUT_DIR}/${PATIENT}_RG.bam"
DEDUP_BAM="${OUTPUT_DIR}/${PATIENT}_dedup.bam"
SPLIT_BAM="${OUTPUT_DIR}/${PATIENT}_split.bam"

# Define Final GVCF Output
GVCF_OUTPUT="${GVCF_DIR}/${PATIENT}_variants.g.vcf.gz"

cd $OUTPUT_DIR

echo "========================================="
echo "STARTING GVCF PIPELINE FOR: $PATIENT"
echo "========================================="

# --- 2. PREPROCESSING PIPELINE ---

# A. Add Read Groups
gatk AddOrReplaceReadGroups \
    -I $STAR_BAM \
    -O $RG_BAM \
    -SO coordinate \
    -RGID $PATIENT -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $PATIENT

# B. Mark Duplicates
gatk MarkDuplicates \
    -I $RG_BAM \
    -O $DEDUP_BAM \
    -M ${PATIENT}_dup_metrics.txt

# C. SplitNCigarReads (Splits reads at splice junctions for RNA-seq)
gatk SplitNCigarReads \
    -R $REF \
    -I $DEDUP_BAM \
    -O $SPLIT_BAM

# --- 3. VARIANT CALLING (GVCF MODE) ---
gatk HaplotypeCaller \
    -R $REF \
    -I $SPLIT_BAM \
    -O $GVCF_OUTPUT \
    -ERC GVCF \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20

# --- 4. CLEANUP ---
# Remove intermediate BAM files to conserve storage space
rm $RG_BAM $DEDUP_BAM $SPLIT_BAM

echo "Pipeline complete for $PATIENT. GVCF successfully saved."
