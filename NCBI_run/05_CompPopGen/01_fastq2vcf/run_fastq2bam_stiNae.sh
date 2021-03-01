#!/bin/bash
#SBATCH -J sm_stiNae
#SBATCH -o out_stiNae
#SBATCH -e err_stiNae
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_fastq2bam_stiNae --profile ./profiles/slurm
