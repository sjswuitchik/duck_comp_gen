#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o logs/slurm-%j.out
#SBATCH -e logs/slurm-%j.err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00-24:00:00
#SBATCH --mem=8000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

source activate raxml

for file in trimmed/subset/*.fa;
do
  raxml-ng --redo --msa $file --model HKY+G4 --prefix $file
done
