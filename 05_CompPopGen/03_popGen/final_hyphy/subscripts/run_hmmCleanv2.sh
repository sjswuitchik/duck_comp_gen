#!/bin/bash
#SBATCH -J hmmCl
#SBATCH -o out_hmv2
#SBATCH -e err_hmv2
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/final_hyphy/aligned

for file in *.filtered;
do
if [ -f ${file}.filtered_hmm.fasta]; then
  continue;
fi
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl $file
done
