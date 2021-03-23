#!/bin/bash

# run from /n/holyscratch01/informatics/swuitchik/snakemake/hetAtr_run

# ./vcfHeaderClean.sh

module load htslib/1.5-fasrc02 bcftools/1.5-fasrc02

zcat hetAtr.vcf.gz | \
sed -e '/^##FILTER=<ID=LowQual,Description=/d' \
-e '/^##FILTER=<ID=GATK_default,Description=/d' \
-e '/^##FILTER=<ID=FS_SOR_filter,Description=/s/.*/##FS_SOR_filter,Description="filter SNPs for strand bias with Phred-scaled p-value for Fishers exact test above 60 and Symmetric Odds Ratio above 3; or indels; or if mixed, a Phred-scaled p-value above 200 and a Symmetric Odds Ratio above 10"/' \
-e '/^##FILTER=<ID=MQ_filter,Description=/s/.*/##MQ_filter,Description="filter SNPs with RMS mapping quality less than 40 and Z-score for Wilcoxon rank sum test for read mapping quality less than -12.5"/' \
-e '/^##FILTER=<ID=RPRS_filter,Description=/s/.*/##RPRS_filter,Description="filter SNPs with Z-Score for Wilcoxon rank sum test for read position bias less than -8; or indels; or if mixed, a Z-Score for Wilcoxon rank sum test less than -20"/' | \
bgzip -c > hetAtr.clean.vcf.gz

bcftools index -t hetAtr.clean.vcf.gz

zcat stiNae.vcf.gz | \
sed -e '/^##FILTER=<ID=LowQual,Description=/d' \
-e '/^##FILTER=<ID=GATK_default,Description=/d' \
-e '/^##FILTER=<ID=FS_SOR_filter,Description=/s/.*/##FS_SOR_filter,Description="filter SNPs for strand bias with Phred-scaled p-value for Fishers exact test above 60 and Symmetric Odds Ratio above 3; or indels; or if mixed, a Phred-scaled p-value above 200 and a Symmetric Odds Ratio above 10"/' \
-e '/^##FILTER=<ID=MQ_filter,Description=/s/.*/##MQ_filter,Description="filter SNPs with RMS mapping quality less than 40 and Z-score for Wilcoxon rank sum test for read mapping quality less than -12.5"/' \
-e '/^##FILTER=<ID=RPRS_filter,Description=/s/.*/##RPRS_filter,Description="filter SNPs with Z-Score for Wilcoxon rank sum test for read position bias less than -8; or indels; or if mixed, a Z-Score for Wilcoxon rank sum test less than -20"/' | \
bgzip -c > stiNae.clean.vcf.gz

bcftools index -t stiNae.clean.vcf.gz
