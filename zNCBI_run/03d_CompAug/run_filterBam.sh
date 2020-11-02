#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p shared
#SBATCH -e filterBam_%A.e
#SBATCH -o filterBam_%A.o
#SBATCH --mem=50000  
#SBATCH -t 23:00:00

singularity shell --cleanenv augustus-2020-05-27-1b69b25ed001.sif
filterBam --in oxyJam.input.sort.bam --out oxyJam.filter.bam --uniq 
