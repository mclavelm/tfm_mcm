# tfm_mcm
**Bioinformatics pipeline for analyzing AmpliSeq data from in vitro adipocytes derived from mesenchymal stem cells and in vivo adipose tissue.** Includes quality control, trimming, alignment, and differential expression analysis. Part of a Master's Thesis. 
## Overview
This repository contains scripts for processing sequencing data:
1. **Trimming** (`trimming.sh`)
2. **Alignment** (`mapping_genome.sh`, `mapping_trans.sh`)
3. **Differential expression analysis with salmon** (`dea_salmon.qmd`)
4. **Clinical data-based analysis** (`dea_datos_clinicos.qmd`)
5. **Co-expression modules analysis** (`modulos_coexpresion.qmd`)  
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
This .qmd file performs differential expression analysis on the quantification results obtained using the `salmon.sh` script.      
**Dependencies**   
-[`edgeR`] (https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
-[`ggplot2`] (https://ggplot2.tidyverse.org/)  
-[`clusterProfiler`] (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)  
-[`org.Hs.eg.db`] (https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)   

## **4. Clinical data-based analysis: `dea_datos_clinicos.qmd`**
This .qmd file integrates clinical data from patients with the expression data obtained from the previous steps.  
It attempts to correlate clinical variables with the gene expression profiles to explore potential patterns.  
It serves as a basis for further refinement in incorporating clinical data into differential expression analysis.  

## **5. Co-expression modules analysis: `modulos_coexpresion.qmd`**  
This .qmd file focuses on identifying gene co-expression modules using the WGCNA  (Weighted Gene Co-expression Network Analysis) approach.   
**Dependencies**   
-[`WGCNA`] (https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html)  
-[`edgeR`] (https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
-[`ggplot2`] (https://ggplot2.tidyverse.org/)  
-[`clusterProfiler`] (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)  
-[`org.Hs.eg.db`] (https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)  



