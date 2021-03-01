#!/bin/bash
#SBATCH -J sm_netAur
#SBATCH -o out_netAur
#SBATCH -e err_netAur
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_fastq2bam_netAur --profile ./profiles/slurm
