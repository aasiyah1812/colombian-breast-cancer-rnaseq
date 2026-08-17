# Unravelling Ancestry-Specific Genomic Differences in Breast Cancer Using RNA-Seq

Welcome to the code repository for my Master's thesis in Molecular Biosciences at Goethe University Frankfurt. 

## 📌 Project Overview
This repository contains the complete bioinformatics pipeline and statistical R scripts used to investigate the biological impact of genetic ancestry and its influence on neoadjuvant chemotherapy response in a highly admixed Colombian breast cancer cohort (N=58). 

The computational workflow dynamically infers genetic ancestry directly from active tumor RNA sequencing data and contrasts it against static germline DNA baselines to uncover transcriptomic masking effects. It also includes all downstream differential gene expression (DESeq2) and functional pathway analyses.

## 🧬 Data Provenance and Acknowledgements
The raw RNA sequencing data and foundational clinical metadata processed by these scripts were generously provided by Dr. Michelle Guevara-Nieto and colleagues. The dataset represents a unique, highly admixed Colombian cohort and is publicly available in the NCBI Gene Expression Omnibus (GEO) under accession number **GSE280902**. 

For the comprehensive clinical and epidemiological framework of this patient cohort, please refer to their foundational publication: *Guevara-Nieto et al. (2025)*.

## 📂 Repository Structure & Thesis Mapping
The scripts in this repository are numbered sequentially with an "A" prefix (e.g., A.1, A.2) to perfectly match the Appendix references cited throughout the Materials and Methods section of my thesis manuscript. 

* **A1 - A3 (Pre-processing & Alignment):** Bash scripts for merging raw FASTQ files, running FastQC/MultiQC, and executing STAR alignment on the Fuchs High-Performance Computing cluster.
* **A4 - A9 (Variant Calling & Population Genomics):** The GATK pipeline (HaplotypeCaller in GVCF mode), 1000 Genomes Project intersection, dataset merging, and ADMIXTURE ancestry estimation.
* **A10 - A16 (Visualization & Divergence Statistics):** R scripts for PCA plotting, evaluating ADMIXTURE cross-validation errors, and quantifying the transcriptomic vs. germline ancestry divergence (Wilcoxon signed-rank tests).
* **A17 - A26 (Transcriptomics & Clinical Analysis):** R scripts for cohort stratification, DESeq2 differential expression models, clinical correlation statistics (Monte Carlo simulated Fisher's Exact tests), and Gene Ontology/GSEA pathway enrichment.

## 💻 Software Prerequisites
To ensure computational reproducibility, the following major software suites were utilized in this pipeline:
* **Bash Environment:** STAR aligner, GATK4, BCFtools, PLINK, ADMIXTURE, FastQC, MultiQC.
* **R Environment (v4.3.2):** tidyverse, DESeq2, clusterProfiler, ggplot2, pheatmap, ggrepel.

## 📬 Contact
If you have any questions regarding the pipeline, specific parameter choices, or the thesis research itself, please feel free to reach out.
