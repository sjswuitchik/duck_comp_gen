#!/bin/sh
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8000
#SBATCH -t 1-00:00:00

module load bedtools2/2.26.0-fasrc01

bedtools unionbedg -empty -g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/hetAtr/genome/hetAtr.fa -i hetAtr_ind*.statscov.bg.gz > hetAtr_union.bg

gzip hetAtr_union.bg
