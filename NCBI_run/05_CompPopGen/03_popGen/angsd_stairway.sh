# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_stiNae_qc

source activate vcfqc
# get Stairway
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
unzip stairway_plot_v2.1.1.zip
rm stairway_plot_v2.1.1.zip
mv stairway-plot-v2/ stairway/

# ANGSD
sbatch run_angsd.sh
angsd -b bams -GL 2 -out hetAtr -doMajorMinor -doMaf -remove_bads 1 -checkBamHeaders 0
angsd -b bams -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -anc hetAtr.ncbi.fasta -ref hetAtr.ncbi.fasta -out hetAtr
realSFS hetAtr.saf.idx > hetAtr.sfs



cp ../hetAtr.filtered.vcf .
cp ../hetAtr_indvs .
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/genomes/hetAtr.ncbi.fasta .
samtools faidx hetAtr.ncbi.fasta -o hetAtr.ncbi.fasta.fai
vcftools --vcf hetAtr.filtered.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --out hetAtr.stair --recode --recode-INFO-all
angsd -vcf-gl hetAtr.stair.recode.vcf -doSaf 1 -out hetAtr -anc hetAtr.ncbi.fasta 
realSFS hetAtr.saf.idx -fold 1 -P 16 > hetAtr.folded



#!/bin/sh
#SBATCH -p bigmem
#SBATCH --mem=190g
#SBATCH -c 16 
#SBATCH --time=24:00:00
#SBATCH --mail-user=nweeks@g.harvard.edu
#SBATCH --mail-type=END,FAIL
source activate vcfqc
realSFS hetAtr.saf.idx -fold 1 -P ${SLURM_JOB_CPUS_PER_NODE} > hetAtr.sfs
