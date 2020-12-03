#!/bin/sh
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem=64000
#SBATCH -p shared
#SBATCH -e starprep.e
#SBATCH -o starprep.o

module load Anaconda/5.0.1-fasrc01
source activate compAug

STAR --runMode genomeGenerate --runThreadN ${SLURM_JOB_CPUS_PER_NODE} --genomeDir genome/ --genomeFastaFiles genome/oxyJam.fa --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 13
