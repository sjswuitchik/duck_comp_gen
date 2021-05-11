#!/bin/bash
#SBATCH -J sm_stiNae
#SBATCH -o out_stiNae
#SBATCH -e err_stiNae
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_intervals_stiNae --profile ./profiles/slurm

