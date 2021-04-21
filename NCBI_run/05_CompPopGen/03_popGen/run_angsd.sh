#!/bin/sh
#SBATCH -p bigmem
#SBATCH --mem=190g
#SBATCH -c 16 
#SBATCH --time=24:00:00
source activate vcfqc

angsd -b bams -GL 2 -out hetAtr -doMajorMinor 1 -doMaf 1
angsd -b bams -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -anc hetAtr.ncbi.fasta -ref hetAtr.ncbi.fasta -out hetAtr
realSFS hetAtr.saf.idx -fold 1 -P ${SLURM_JOB_CPUS_PER_NODE} > hetAtr.sfs
