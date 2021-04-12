#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 00-24:00:00
#SBATCH --mem=8000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

source activate raxml

for file in trimmed/*.fa;
do
  raxml-ng --msa $file --model HKY+G4 --prefix $file --seed 2 --redo
done

cat trimmed/*.bestTree > trimmed/final.tree

