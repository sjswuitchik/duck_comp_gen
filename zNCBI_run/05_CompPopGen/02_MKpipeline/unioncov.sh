#!/bin/sh
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8000
#SBATCH -t 1-00:00:00

module load bedtools2/2.26.0-fasrc01

bedtools unionbedg -header -empty-g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/hetAtr/genome/hetAtr.fa -i hetAtr_ind01.statscov.bg hetAtr_ind02.statscov.bg hetAtr_ind03.statscov.bg hetAtr_ind04.statscov.bg hetAtr_ind05.statscov.bg hetAtr_ind06.statscov.bg hetAtr_ind07.statscov.bg hetAtr_ind08.statscov.bg hetAtr_ind09.statscov.bg hetAtr_ind10.statscov.bg hetAtr_ind11.statscov.bg hetAtr_ind12.statscov.bg hetAtr_ind13.statscov.bg hetAtr_ind14.statscov.bg > hetAtr_union.bg

gzip hetAtr_union.bg
