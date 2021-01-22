#!/usr/bin/bash

# VCF filtering
Rscript --vanilla missingness.R $INSHORT'_missing_data.txt' $OUTSHORT'_missing_data.txt'

vcftools --gzvcf $INSHORT.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --remove ingroup.remove.indv --max-missing 0.5 --recode --recode-INFO-all --out $INSHORT.filter
vcftools --gzvcf $OUTSHORT.vcf.gz --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf 0 --remove outgroup.remove.indv --max-missing 0.5 --recode --recode-INFO-all --out $OUTSHORT.filter

# create callable sites for in and outgroup
bedtools intersect -a $INSHORT'_coverage_sites_clean_merged.bed' -b $OUTSHORT'_coverage_sites_clean_merged.bed' > callable.bed
# intersect with in and out group to filter for callable sites only
bedtools intersect -a $INSHORT.filter.recode.vcf -b callable.bed -header > $INSHORT.call.vcf 
bedtools intersect -a $OUTSHORT.filter.recode.vcf -b callable.bed -header > $OUTSHORT.call.vcf

# annotates ingroup VCF
java -jar $PATHS/snpEff.jar $INSHORT $PATHW/$INSHORT.call.vcf > $PATHW/$INSHORT.ann.vcf
# annotate outgroup VCF
java -jar $PATHS/snpEff.jar $INSHORT $PATHW/$OUTSHORT.call.vcf > $PATHW/$OUTSHORT.ann.vcf

#parse out variants from annotated VCF and output an annotated BED file 
python annot_parser.py $INSHORT.ann.vcf $INSHORT.ann.bed -key missense_variant -key synonymous_variant
python annot_parser.py $OUTSHORT.ann.vcf $OUTSHORT.ann.bed -key missense_variant -key synonymous_variant

# pull out only CDS regions
column -s, -t < genes.gff | awk '$3 == "CDS"' > onlyCDS.gff
# convert to bed
awk -f gff2bed.awk onlyCDS.gff > onlyCDS.bed
# pull out gene names
cat onlyCDS.bed | python genenames.py > onlyCDS.genes.bed

# associate gene regions with annotations
bedtools intersect -a $INSHORT.ann.bed -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > $INSHORT.final.bed
bedtools intersect -a $OUTSHORT.ann.bed -b onlyCDS.genes.bed -wb | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > $OUTSHORT.final.bed

# associate callable sites with genes
bedtools intersect -a callable.bed -b onlyCDS.genes.bed -wb | cut -f1,2,3,7 | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct > callable.cds.bed 
# output site missingness information for use in SnIPRE
vcftools --vcf $INSHORT.ann.vcf --missing-site --out $INSHORT.missingness
vcftools --vcf $OUTSHORT.ann.vcf --missing-site --out $OUTSHORT.missingness

# prep data for SnIPRE, run SnIPRE, run MK test, and calculate direction of selection (DOS)
Rscript --slave --vanilla mktest.R $INSHORT'.final.bed' $OUTSHORT'.final.bed' $INSHORT'.missingness.lmiss' $OUTSHORT'.missingness.lmiss' >std.Rout 2>std.Rerr

