#!/bin/bash
# --- SBATCH Directives ---
#SBATCH --job-name=Download_1kG_Phased
#SBATCH --partition=fuchs
#SBATCH --output=/path/to/project_directory/logs/Download_1kG_%j.out
#SBATCH --error=/path/to/project_directory/logs/Download_1kG_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

# --- 1. ENVIRONMENT & PATHS ---
echo "Starting 1000 Genomes High Coverage Phased Data Download..."

# Set up directory for the global reference panel
BASE_DIR="/path/to/project_directory"
REF_DIR="${BASE_DIR}/1000G_phased_GRCh38"

mkdir -p ${REF_DIR}
cd ${REF_DIR}

# --- 2. AUTOMATED FILE RETRIEVAL ---
# Loop through autosomes 1 to 22 to download phased VCFs and index files
for chr in {1..22}
do
    echo "Downloading GRCh38 phased data for chr ${chr}..."
    
    # Download the compressed VCF file (-c flag allows resuming interrupted downloads)
    wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
    
    # Download the corresponding tabix index file
    wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
done

echo "Download of all autosomal reference data complete."
