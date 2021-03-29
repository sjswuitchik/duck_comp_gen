#!/bin/bash
#SBATCH -J muscleAlign
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=8000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

conda activate align

for file in fastas/*.fa;
do
  muscle -in $file -quiet -out $file.afa
done
