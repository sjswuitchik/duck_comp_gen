#!/bin/bash
#SBATCH -J RAxMLsub
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00-03:00:00
#SBATCH --mem=4000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

for file in fastas/muscle/subset/*.afa;
do
  raxml-ng --msa $file --model HKY+G4 --prefix T_$file --seed 2
done
