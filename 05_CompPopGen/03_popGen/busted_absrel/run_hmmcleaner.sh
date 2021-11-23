#!/bin/bash
#SBATCH -J hmmCl
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=15000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft

for file in *.mafft;
do
if [ -f $file_hmm.fasta]; then
continue;
fi
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl $file
done
