#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 8
#SBATCH -t 06-23:00:00
#SBATCH --mem=12000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/trimmed/finished_fastas/

source activate raxml

for file in *.fa;
do
  raxml-ng --msa $file --model HKY+G4 --prefix $file --all --seed 2 --redo --bs-metric fbp
done

cat *.bestTree > final.all.tree

cat *.support > final.support.tree
