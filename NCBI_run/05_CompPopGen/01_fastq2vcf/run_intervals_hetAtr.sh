#!/bin/bash
#SBATCH -J sm_hetAtr
#SBATCH -o out_hetAtr
#SBATCH -e err_hetAtr
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000


source activate snakemake
snakemake --snakefile Snakefile_intervals_hetAtr --profile ./profiles/slurm

