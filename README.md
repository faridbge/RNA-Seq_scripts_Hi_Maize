
# Code Pipeline for RNA-Seq Mapping and Read Quantification in Hi-A Corn

Welcome to this repository! This collection showcases the codes and scripts used for the Hi-A Corn tanscriptomics project. The repository includes scripts for RNA-Seq mapping and read quantification. All these analyses were conducted in the Texas A&M High Performance Research Computing (HPRC) cluster.

Here's what you'll find:
- **README file: describe about the page**
- **Script for RNA mapping:**
    - *RNA Mapping Pipeline.md*
    - *fastqc_terrra.slurm* (Quality check of raw reads)
    - *trimmomatic_loop_terra.slurm* (Trimming reads)
    - *trimmomatic_terra.slurm* (Trimming reads, alternative script)
    - ref_genome_index_star_terra.slurm (Indexing reference genome)
    - mapping_reads_indexed_ref_genome.slurm (Mapping reads to reference genome)
      
- **Script for reads counting**:
     
    - *Counting reads*
        - script_salmon_quantifying.slurm
        - script_salmon_indexing.slurm
        - script_counting_STAR_summary   
    - *Analysis:*
        - script_processing_STAR_output_files.R
          
            
## Languages Used

- Bash
- R

