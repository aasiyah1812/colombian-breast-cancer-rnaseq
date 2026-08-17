#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=STAR_Align_Array
#SBATCH --partition=fuchs
#SBATCH --array=1-10
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --output=logs/STAR_align_%A_%a.out
#SBATCH --error=logs/STAR_align_%A_%a.err

# --- 1. Define Patient ID from Array Index ---
# Extracts the specific patient ID based on the array task number
PATIENT_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" patient_list.txt)

if [ -z "$PATIENT_ID" ]; then
    echo "Error: No patient ID provided by array list."
    exit 1
fi

echo "Processing Patient: $PATIENT_ID"

# --- 2. SYSTEM PATHS & REFERENCES ---
REF="/path/to/reference/GRCh38.p14.genome.fa"
GTF="/path/to/annotation/gencode.v48.annotation.gtf"
STAR_INDEX_DIR="/path/to/reference/STAR_Index"
STAR_EXEC="/path/to/software/STAR-2.7.11b/bin/Linux_x86_64/STAR"
THREADS=${SLURM_CPUS_PER_TASK}

# --- 3. DIRECTORIES ---
BASE_DIR="/path/to/project_directory"
RAW_DATA_DIR="${BASE_DIR}/RAW_FASTQ_DATA"
OUTPUT_DIR="${BASE_DIR}/${PATIENT_ID}_analysis"

mkdir -p ${BASE_DIR}/logs
mkdir -p ${OUTPUT_DIR}

# --- 4. ACTIVATE VIRTUAL ENVIRONMENT ---
source /path/to/virtualenv/bin/activate

# --- 5. PREPARE INPUT FILES ---
R1_TARGET="${OUTPUT_DIR}/${PATIENT_ID}_R1_merged.fastq.gz"
R2_TARGET="${OUTPUT_DIR}/${PATIENT_ID}_R2_merged.fastq.gz"

echo "--- STEP 1: Preparing Input Files ---"

RAW_R1="${RAW_DATA_DIR}/${PATIENT_ID}_R1_merged.fastq.gz"
RAW_R2="${RAW_DATA_DIR}/${PATIENT_ID}_R2_merged.fastq.gz"

# Check if merged file exists; if not, concatenate raw lanes dynamically
if [ -f "$RAW_R1" ]; then
    echo "Found pre-merged files. Creating symbolic links..."
    ln -sf $RAW_R1 $R1_TARGET
    ln -sf $RAW_R2 $R2_TARGET
else
    echo "Pre-merged file not found. Concatenating available lanes..."
    cat ${RAW_DATA_DIR}/${PATIENT_ID}_*_R1_*.fastq.gz > ${R1_TARGET}
    cat ${RAW_DATA_DIR}/${PATIENT_ID}_*_R2_*.fastq.gz > ${R2_TARGET}
fi

# --- 6. STAR ALIGNMENT & QUANTIFICATION ---
PREFIX="${PATIENT_ID}_star"

echo "--- STEP 2: STAR Alignment ---"
${STAR_EXEC} --runThreadN ${THREADS} \
     --genomeDir ${STAR_INDEX_DIR} \
     --readFilesIn ${R1_TARGET} ${R2_TARGET} \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPUT_DIR}/${PREFIX}_ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile ${GTF} \
     --outFilterMultimapNmax 20

echo "STAR Alignment Complete for ${PATIENT_ID}."
