#!/bin/bash
#SBATCH -J gatkUpdate
#SBATCH -o out_%j
#SBATCH -e err_%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=8000

# for duck filtering update
# sbatch run_gatkUpdate.sh Combined_hardFiltered spp_name

set -o errexit

source activate gatk

picard SortVcf -I $1.vcf.gz -O $1_sorted.vcf.gz

gatk IndexFeatureFile -I $1_sorted.vcf.gz

gatk VariantFiltration -R /n/holyscratch01/informatics/swuitchik/snakemake/$2_run/genome/$2.fa -V $1_sorted.vcf.gz --output $1_updatedFilter.vcf.gz --filter-name "RPRS_filter" --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)" --filter-name "FS_SOR_filter" --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" --filter-name "MQ_filter" --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))"
