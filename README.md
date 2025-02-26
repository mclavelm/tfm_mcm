# tfm_mcm
**Bioinformatics pipeline for analyzing AmpliSeq data from mesenchymal stem cells and adipose tissue.** Includes quality control, trimming, alignment, and differential expression analysis. Part of a Master's Thesis. 
## Overview
This repository contains scripts for processing sequencing data:
1. **Trimming** (`trimming.sh`)
2. **Alignment** (`mapping_genome.sh`, `mapping_trans.sh`)
3. **Alignment statistics** (`alignment_stats.sh`)

##**1. Trimming: `trimming.sh`**
This script performs quality trimming on raw FASTQ files using fastp.
**Dependencies**
[`fastp`] (https://github.com/OpenGene/fastp)
Run the script inside the folder containing raw FASTQ files.
