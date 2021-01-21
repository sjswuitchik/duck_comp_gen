#!/bin/sh
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8000
#SBATCH -t 1-00:00:00

module load bedtools2/2.26.0-fasrc01

bedtools unionbedg -empty -g hetAtr.chrom.sizes -i hetAtr_ind*.statscov.bg > hetAtr_union.bg

gzip hetAtr_union.bg
