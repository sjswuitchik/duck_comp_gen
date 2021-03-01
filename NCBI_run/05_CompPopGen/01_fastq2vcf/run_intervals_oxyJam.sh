#!/bin/bash
#SBATCH -J sm_oxyJam
#SBATCH -o out_oxyJam
#SBATCH -e err_oxyJam
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

module load Anaconda3/2020.11
source activate snakemake
snakemake --snakefile Snakefile_intervals_oxyJam --profile ./profiles/slurm
