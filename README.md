
# Code Pipeline for RNA-Seq Mapping and Read Quantification in Hi-A Corn

Welcome to this repository! This collection contains the scripts and workflows used the Hi-A Corn Transcriptomics Project. These analyses were conducted on the Texas A&M High Performance Research Computing (HPRC) cluster.

## Overview
This repository provides a step-by-step workflow for processing RNA-Seq data, from raw read quality control to reads quantification
## Directory Structure
**01_mapping:** This section covers the alignment of raw reads to the reference genome using STAR.
- *fastqc_terrra.slurm*: Quality check of raw reads
- *trimmomatic_loop_terra.slurm*: Trimming of adapter sequences
- *trimmomatic_terra.slurm* (Trimming reads, alternative script)
- *ref_genome_index_star_terra.slurm*: Generation of the STAR reference genome index
- *mapping_reads_indexed_ref_genome.slurm*: Mapping reads to the indexed reference genome
    
**02_quantification**: This folder contains scripts for reads counting.     
    - *script_counting_STAR_summary*
    - *script_salmon_indexing.slurm*
    - *script_salmon_quantifying.slurm*
    - *script_processing_STAR_output_files.R*              
            
## Languages Used
- Bash
- R

