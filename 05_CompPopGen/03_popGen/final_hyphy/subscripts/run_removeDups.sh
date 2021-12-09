#!/bin/bash
#SBATCH -J removeDups
#SBATCH -o out_dups
#SBATCH -e err_dups
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/

source activate align

while IFS= read -r file
do
  hyphy remove-duplicates.bf --msa aligned/${file}.clean.fa -- tree gene_trees/${file}_tree.txt
done < "clean_ogs.tsv" 
