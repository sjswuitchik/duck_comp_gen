#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o logs/slurm-%j.out
#SBATCH -e logs/slurm-%j.err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=12000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

source activate raxml

for file in trimmed/*.fa;
do
  raxml-ng --msa $file --model HKY+G4 --prefix $file --seed 2 --redo
done

cat trimmed/*.bestTree > trimmed/final.tree

