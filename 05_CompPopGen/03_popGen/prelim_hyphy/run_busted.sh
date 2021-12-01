#!/bin/bash
#SBATCH -J busted
#SBATCH -o out_bus
#SBATCH -e err_bus
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/clean_align/all_spp/

source activate align

for file in *.filtered;
do
if [ -f $file.BUSTED.json ]; then
continue;
fi
hyphy busted --alignment $file --tree gallo.newick
done
