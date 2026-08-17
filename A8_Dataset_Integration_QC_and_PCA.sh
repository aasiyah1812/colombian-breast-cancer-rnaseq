#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=Merge_QC_PCA
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/Merge_QC_%j.out
#SBATCH --error=/path/to/project_directory/logs/Merge_QC_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# --- 1. ENVIRONMENT & PATHS ---
# Activate virtual environment
source /path/to/virtualenv/bin/activate

BASE_DIR="/path/to/project_directory"
OUT_DIR="${BASE_DIR}/05_Pipeline_Results"
REF_1000G="${BASE_DIR}/02_Reference_Data/True_1000G_Reference.vcf.gz"

PLINK="/path/to/software/plink"

cd ${OUT_DIR}

# --- 2. DATASET CONVERSION ---
echo "Step 1: Converting VCFs to PLINK binary format..."
${PLINK} --vcf cohort_58_intersected.vcf.gz --make-bed --out plink_cohort_58 \
    --set-missing-var-ids '@:#:$1:$2' --double-id

# Standardize nomenclature: Remove "chr" prefix to identically match the 1000G database
sed -i 's/chr//g' plink_cohort_58.bim

# Convert 1000G reference (restricted to SNPs)
${PLINK} --vcf ${REF_1000G} --snps-only --make-bed --out plink_ref_1000G \
    --set-missing-var-ids '@:#:$1:$2' --double-id

# --- 3. VARIANT HARMONIZATION ---
echo "Step 2: Identifying and extracting consensus SNPs..."
# Identify overlapping variants between the cohort and the reference panel
awk '{print $2}' plink_cohort_58.bim > overlapping_snps.txt
${PLINK} --bfile plink_ref_1000G --extract overlapping_snps.txt --make-bed --out plink_ref_1000G_subset

# Extract strictly overlapping markers to ensure exact dataset congruence prior to merging
awk '{print $2}' plink_ref_1000G_subset.bim > true_overlapping_snps.txt
${PLINK} --bfile plink_cohort_58 --extract true_overlapping_snps.txt --make-bed --out plink_cohort_58_subset

# --- 4. DATASET MERGING ---
echo "Step 3: Merging the cohort with the reference panel..."
rm -f merged_raw-merge.missnp

# Execute primary merge attempt
${PLINK} --bfile plink_cohort_58_subset --bmerge plink_ref_1000G_subset --make-bed --out merged_raw

# Automated handling for allelic/strand mismatches
if [ -f merged_raw-merge.missnp ]; then
    echo "Conflicting strand SNPs detected. Excluding mismatched variants..."
    ${PLINK} --bfile plink_cohort_58_subset --exclude merged_raw-merge.missnp --make-bed --out plink_cohort_58_clean
    echo "Retrying dataset merge..."
    ${PLINK} --bfile plink_cohort_58_clean --bmerge plink_ref_1000G_subset --make-bed --out merged_raw
fi

# --- 5. POPULATION GENETICS QUALITY CONTROL ---
echo "Step 4: Applying quality control filters..."
# Filter for MAF > 0.05, variant missingness < 10%, individual missingness < 10%
${PLINK} --bfile merged_raw --maf 0.05 --geno 0.1 --mind 0.1 --make-bed --out merged_qc

# --- 6. LINKAGE DISEQUILIBRIUM PRUNING ---
echo "Step 5: Performing LD pruning for independent markers..."
${PLINK} --bfile merged_qc --indep-pairwise 50 10 0.1 --out prune_data
${PLINK} --bfile merged_qc --extract prune_data.prune.in --make-bed --out FINAL_Ancestry_Data

# --- 7. PRINCIPAL COMPONENT ANALYSIS ---
echo "Step 6: Computing Principal Component Analysis..."
${PLINK} --bfile FINAL_Ancestry_Data --pca --out PCA_Results

echo "Integration and Population Genetics Pipeline Complete."
