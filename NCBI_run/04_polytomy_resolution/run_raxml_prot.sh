#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o logs/out_%j
#SBATCH -e logs/err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=12000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/raxml_prot_code

conda activate raxml

for file in *.filtered;
do
  raxml-ng --msa $file --model HKY+G4 --prefix $file --seed 2 --redo
done

cat *.bestTree > final.tree
