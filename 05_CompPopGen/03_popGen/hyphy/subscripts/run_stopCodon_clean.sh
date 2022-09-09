#!/bin/bash
#SBATCH -J stopCodon
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/aligned

for file in *_hmm.fasta;
do
  sed 's/tag/\-\-\-/g' $file > $file.filt1
  sed 's/taa/\-\-\-/g' $file.filt1 > $file.filt2
  sed 's/tga/\-\-\-/g' $file.filt2 > $file.filtered
  rm $file.filt1
  rm $file.filt2
done
