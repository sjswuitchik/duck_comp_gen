#!/bin/bash
#SBATCH -J busted
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned

source activate align

for file in *.fa_hmm.fasta;
do
hyphy busted --alignment $file --tree gallo.newick
done
