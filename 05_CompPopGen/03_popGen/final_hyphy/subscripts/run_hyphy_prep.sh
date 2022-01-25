#!/bin/bash
#SBATCH -J hyphyPrep
#SBATCH -o out
#SBATCH -e err
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -t 03-00:00
#SBATCH --mem=9000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/

conda activate align

while IFS= read -r file
do
  hyphy hyphy-analyses/codon-msa/pre-msa.bf --input og_fastas/${file}_nuc.fa
  mafft ${file}_nuc.fa_protein.fas > og_fastas/${file}_protein.msa
  hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa og_fastas/${file}_protein.msa --nucleotide-sequences og_fastas/${file}_nuc.fa_nuc.fas --output og_fastas/${file}_codon.msa --compress No
done < "clean_ogs"
