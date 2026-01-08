
# Code Pipeline for RNA-Seq Mapping and Read Quantification in Hi-A Corn

Welcome to this repository! This collection showcases the codes and scripts used for the Hi-A Corn tanscriptomics project. The repository includes scripts for RNA-Seq mapping and read quantification. All these analyses were conducted in the Texas A&M High Performance Research Computing (HPRC) cluster.

Here's what you'll find:
- **README file**
- **Script documentation:**
    - *RNA Mapping Pipeline.md*
      
- **Scripts:**
    - *Preprocessing:*
        - fastqc_terrra.slurm (Quality check of raw reads)
        - trimmomatic_loop_terra.slurm (Trimming reads)
        - trimmomatic_terra.slurm (Trimming reads, alternative script)
    - *Alignment:*
        - ref_genome_index_star_terra.slurm (Indexing reference genome)
        - mapping_reads_indexed_ref_genome.slurm (Mapping reads to reference genome)          
    - *Counting reads*
        - script_salmon_quantifying.slurm
        - script_salmon_indexing.slurm
        - script_counting_STAR_summary   
    - *Analysis:*
        - script_processing_STAR_output_files.R
        - script_processing_TPMCalculator_output_files.R
            
## Languages Used

- Bash
- R

