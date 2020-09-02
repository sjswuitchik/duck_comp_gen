#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module load Anaconda/5.0.1-fasrc01
source activate snakemake
snakemake --snakefile Snakefile_bam2vcf_gatk --profile ./profiles/slurm

