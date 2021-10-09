#!/bin/bash
#SBATCH -J absrel
#SBATCH -o out_abs
#SBATCH -e err_abs
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft/clean_align/all_spp

source activate align

for file in *.filtered;
do
if [ -f $file.ABSREL.json ]; then
continue;
fi
hyphy busted --alignment $file --tree gallo.newick
done
