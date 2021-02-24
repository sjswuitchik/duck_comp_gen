#!/bin/bash
#SBATCH -J sm_hetAtr
#SBATCH -o out_hetAtr
#SBATCH -e err_hetAtr
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module load Anaconda3/2020.11
source activate snakemake
snakemake --snakefile Snakefile_fastq2bam_hetAtr --profile ./profiles/slurm
