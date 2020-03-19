#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -p serial_requeue
#SBATCH --mem=2000
#SBATCH --array=0-360

# run from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/alignments/aligned

module load perl/5.26.1-fasrc01 mafft/7.407-fasrc01 trimal/1.2rev59-fasrc01

printf -v BATCH "%03d" $SLURM_ARRAY_TASK_ID
mkdir -p batch${BATCH}_output

for CNEE in $(cat batch$BATCH);
do
  if [ ! -s batch${BATCH}_output/$CNEE.aligned.fa ]
  then 	
    ginsi ../unaligned/$CNEE.fa | perl -p -e 'if (!/^>/) {s/[Nn]/-/g}' > temp_${BATCH}.fa 
    trimal -in temp_${BATCH}.fa -noallgaps -out batch${BATCH}_output/$CNEE.aligned.fa -fasta
    rm temp_${BATCH}.fa
  fi
done

