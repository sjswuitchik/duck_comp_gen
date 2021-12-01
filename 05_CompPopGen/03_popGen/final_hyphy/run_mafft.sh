#!/bin/bash
#SBATCH -J mafft
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/final_hyphy/

source activate align

for file in *.fa;
do
if [ -f $file.mafft ]; then
  continue;
fi
  mafft --quiet og_fastas/${file} > aligned/${file}.mafft
done
