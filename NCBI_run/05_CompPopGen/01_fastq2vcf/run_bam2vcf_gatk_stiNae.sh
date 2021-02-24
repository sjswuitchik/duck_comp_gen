#!/bin/bash
#SBATCH -J sm_stiNae
#SBATCH -o out_stiNae
#SBATCH -e err_stiNae
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module load Anaconda3/2020.11
source activate snakemake
snakemake --snakefile Snakefile_bam2vcf_gatk_stiNae --profile ./profiles/slurm

