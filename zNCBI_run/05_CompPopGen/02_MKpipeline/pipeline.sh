#!/usr/bin/bash

# VCF concatenation

cd $INSHORT.vcfs 
bcftools concat *.vcf.gz -O v -o $INSHORT.concat.vcf
mv $INSHORT.concat.vcf ..
cd ../$OUTSHORT.vcfs 
bcftools concat *.vcf.gz -O v -o $OUTSHORT.concat.vcf
mv $OUTSHORT.concat.vcf ..
cd ..

# VCF filtering
vcftools --vcf $INSHORT.concat.vcf --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --mac 1 --max-missing 0.75 --recode --recode-INFO-all --out $INSHORT.clean
mv $INSHORT.clean.recode.vcf $INSHORT.clean.vcf

vcftools --vcf $OUTSHORT.concat.vcf --remove-filtered-all --remove-indels --min-alleles 2 --max-alleles 2 --maf 0 --max-missing 0.75 --recode --recode-INFO-all --out $OUTSHORT.clean
mv $OUTSHORT.clean.recode.vcf $OUTSHORT.clean.vcf

Rscript --vanilla missingness.R $INLONG'_all_all_missingness_info.txt' $OUTLONG'_all_all_missingness_info.txt'

# remove individuals with high relative missingness, if any
export ININDV=`cat ingroup.remove.indv | wc -l`
if [ $ININDV -gt 1 ]
then
	vcftools --vcf $INSHORT.clean.vcf --remove-indv ingroup.remove.indv --recode --recode-INFO-all --out $INSHORT.clean2
	mv $INSHORT.clean2.recode.vcf $INSHORT.clean.vcf
else 
	echo "No indv to remove from ingroup"
fi

export OUTINDV=`cat outgroup.remove.indv | wc -l`
if [ $OUTINDV -gt 2 ]
then
	vcftools --vcf $OUTSHORT.clean.vcf --remove-indv outgroup.remove.indv --recode --recode-INFO-all --out $OUTSHORT.clean2
	mv $OUTSHORT.clean2.recode.vcf $OUTSHORT.clean.vcf
else 
	echo "No indv to remove from outgroup"
fi

# create callable sites for in and outgroup
bedtools intersect -a $INLONG'_clean_coverage_sites_merged.bed' -b $OUTLONG'_clean_coverage_sites_merged.bed' > callable.bed
# intersect with in and out group to filter for callable sites only
bedtools intersect -a $INSHORT.clean.vcf -b callable.bed -header > $INSHORT.call.vcf 
bedtools intersect -a $OUTSHORT.clean.vcf -b callable.bed -header > $OUTSHORT.call.vcf

# annotates ingroup VCF
java -jar $PATHS/snpEff.jar $INSHORT $PATHW/$INSHORT.call.vcf > $PATHW/$INSHORT.ann.vcf
# annotate outgroup VCF
java -jar $PATHS/snpEff.jar $INSHORT $PATHW/$OUTSHORT.call.vcf > $PATHW/$OUTSHORT.ann.vcf

#parse out variants from annotated VCF and output an annotated BED file 
python parser_nov.py $INSHORT.ann.vcf $INSHORT.ann.bed -key missense_variant -key synonymous_variant
python parser_nov.py $OUTSHORT.ann.vcf $OUTSHORT.ann.bed -key missense_variant -key synonymous_variant

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

