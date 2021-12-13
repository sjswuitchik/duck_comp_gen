#!/bin/bash
#SBATCH -J hyphyPrep
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 03-00:00
#SBATCH --mem=9000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/og_fastas

conda activate align

for file in *.fa;
do
  hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $file
  muscle -in ${file}_protein.fas -out ${file}_protein.msa
  hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${file}_protein.msa --nucleotide-sequences ${file}_nuc.fas --output ${file}_codon.msa --compress No
done
