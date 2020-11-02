#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p shared
#SBATCH -e star2_%A.e
#SBATCH -o star2_%A.o
#SBATCH -J star
#SBATCH --mem=64000
#SBATCH -t 23:00:00

module load Anaconda/5.0.1-fasrc01
source activate compAug

STAR --runThreadN ${SLURM_JOB_CPUS_PER_NODE} --genomeDir genome --sjdbFileChrStartEnd oxyJamPass1SJ.out.tab --outFileNamePrefix oxyJamPass2 --readFilesManifest mapManifest.tsv
