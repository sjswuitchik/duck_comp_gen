#!/bin/bash
#SBATCH -J sm_netAur
#SBATCH -o out_netAur
#SBATCH -e err_netAur
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module load Anaconda3/2020.11
source activate snakemake
snakemake --snakefile Snakefile_bam2vcf_gatk_netAur --profile ./profiles/slurm

