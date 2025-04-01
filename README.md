# tfm_mcm
**Bioinformatics pipeline for analyzing AmpliSeq data from mesenchymal stem cells and adipose tissue.** Includes quality control, trimming, alignment, and differential expression analysis. Part of a Master's Thesis. 
## Overview
This repository contains scripts for processing sequencing data:
1. **Trimming** (`trimming.sh`)
2. **Alignment** (`mapping_genome.sh`, `mapping_trans.sh`)
3. **Alignment statistics** (`alignment_stats.sh`)

## **1. Trimming: `trimming.sh`**
This script performs quality trimming on raw fastq files using fastp.
**Dependencies**  
-[`fastp`] (https://github.com/OpenGene/fastp)  
Run the script inside the folder containing raw fastq files.

## **2. Quantification: `salmon.sh`**
This script performs the quantification of every fastq file against a previously indexed transcriptome of reference.  
**Dependencies**
-[`salmon`] (https://combine-lab.github.io/salmon/)  
Run the script inside the folder containing the trimmed fastq files. 

## **3. Differential expression analysis and correlated genes: `dea_salmon.qmd`** 

