# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_stiNae_qc
## VariantQC
module load jdk/10.0.1-fasrc01 bcftools/1.5-fasrc02
# get VariantQC tool 
wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.24/DISCVRSeq-1.24.jar
mv DISCVRSeq-1.24.jar DISCVRSeq.jar

# copy over VCF with updated filtering and indv list of hetAtr samples
cp ../hetAtr_run/Combined_hardFiltered_updatedFilter.vcf.gz* ../hetAtr_run/hetAtr_indvs .
# extract indvs and clean up VCFs
bcftools view -O z -S hetAtr_indvs -U -a Combined_hardFiltered_updatedFilter.vcf.gz > hetAtr.filtered.vcf.gz
bcftools view -O z -s stiNae_male -U -a Combined_hardFiltered_updatedFilter.vcf.gz > stiNae.filtered.vcf.gz

# sort
source activate gatk 
picard SortVcf -Xmx8g -I hetAtr.filtered.vcf.gz -O hetAtr.filtered.sorted.vcf.gz
picard SortVcf -Xmx8g -I stiNae.filtered.vcf.gz -O stiNae.filtered.sorted.vcf.gz

# index
gatk IndexFeatureFile -I hetAtr.filtered.sorted.vcf.gz
gatk IndexFeatureFile -I stiNae.filtered.sorted.vcf.gz

# run VariantQC
java -jar DISCVRSeq.jar VariantQC --maxContigs 26704 -R hetAtr.fa -V hetAtr.filtered.sorted.vcf.gz -O hetAtr.filtered.html
java -jar DISCVRSeq.jar VariantQC --maxContigs 26704 -R hetAtr.fa -V stiNae.filtered.sorted.vcf.gz -O stiNae.filtered.html

conda deactivate 

#conda create -n vcfqc -c bioconda plink vcftools bcftools r-base r-tidyverse admixture perl bedtools angsd samtools
source activate vcfqc
# output stats
vcftools --gzvcf hetAtr.filtered.vcf.gz --out hetAtr.rel --relatedness2
# copy callable sites from MK pipeline prep to standardize pi & TajD
cp ../../MKpipeline/hetAtr_coverage_sites_clean_merged.bed .
gunzip hetAtr.filtered.vcf.gz
bedtools intersect -a hetAtr.filtered.vcf -b hetAtr_coverage_sites_clean_merged.bed -header > hetAtr.callable.vcf
vcftools --vcf hetAtr.callable.vcf --out hetAtr.bial.100kb --TajimaD 100000 --min-alleles 2 --max-alleles 2
vcftools --vcf hetAtr.callable.vcf --out hetAtr.pi.bial --window-pi 100000 --min-alleles 2 --max-alleles 2

grep -v '\*' hetAtr.filtered.vcf > hetAtr.filtered.clean.vcf
plink --vcf hetAtr.filtered.clean.vcf --make-bed --geno 0.999 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --indep-pairwise 500 50 0.1 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --make-bed --extract hetAtr.prune.in --out hetAtr.ld_pruned --allow-extra-chr
plink --bfile hetAtr.ld_pruned --ibc --out hetAtr --allow-extra-chr
plink --bfile hetAtr.ld_pruned --pca --out hetAtr --allow-extra-chr

# admixutre
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt
sed 's/\r$//g' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | grep -v "^#" | cut -f1,3,5 > hetAtr_chr_key
awk '{print $3, $2}' hetAtr_chr_key > acckey
# manually editing the chrom names in acckey for now, see if sed or something similar will work later; chr 34 = W, chr 35 = Z, chr 36 = MT
sed -i 's/na/\d0/g' acckey 
./replace_chrs.pl acckey hetAtr.ld_pruned.bim > hetAtr.repl.ld_pruned.bim
mv hetAtr.repl.ld_pruned.bim hetAtr.ld_pruned.bim
for K in {2..5}
do
	admixture --cv hetAtr.ld_pruned.bed $K > hetAtr.${K}.admix.log 2> hetAtr.${K}.admix.err
done
cat hetAtr.*.admix.log | grep CV | perl -pi -e 's/.+=//' | perl -pi -e 's/\): /\t/' > hetAtr.CV

# plot of inbreeding coeficicients 
Rscript hetAtr.plot.r
# plot PCA
Rscript hetAtr.pca.plot
# plot ADMIXTURE results
Rscript hetAtr.admixture.plot.r
