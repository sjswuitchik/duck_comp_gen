#!/bin/bash
#SBATCH -J bedg2bw
#SBATCH -e bedg2bw.err
#SBATCH -o bedg2bw.out
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 02-00:00:00
#SBATCH --mem=10000

# run from /n/holyscratch01/informatics/swuitchik/snakemake/hetAtr_run/fastq2bam_hetAtr/01_mappedReads
# sbatch run_bedg2bw.sh

for file in *.bg.sorted.bg;
do
  sort -k1,1 -k2,2n $file > $file.sorted.bg
  /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/./bedGraphToBigWig $file ../../hetAtr.chrom.sizes $file.bw
done
