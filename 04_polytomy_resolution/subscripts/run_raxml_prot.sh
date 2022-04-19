#!/bin/bash
#SBATCH -J RAxML
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH --exclusive
#SBATCH -t 06-23:00:00
#SBATCH --mem=0

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_coding

source activate raxml

for file in orthos/*.fasta;
do
if [ -f ${file}.raxml.bestTree ]; then
  continue;
fi
  raxml-ng --msa $file --model HKY+G4 --prefix $file --all --seed 2 --redo --bs-metric fbp
done
