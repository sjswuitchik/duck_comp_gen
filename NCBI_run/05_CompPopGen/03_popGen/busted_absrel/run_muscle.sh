#!/bin/bash
#SBATCH -J muscle
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=10000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

conda activate align

for file in *.fa;
do
muscle -in $file -quiet -out $file.afa
done
