# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/vcfs/

module load bcftools/1.5-fasrc02 vcftools/0.1.14-fasrc01 plink/1.90-fasrc01

cp ../shortRead_mapping_variantCalling/gatk/Combined_hardFiltered.vcf .
cp ~/github_repos/duck_comp_gen/NCBI_run/05_CompPopGen/02_MKpipeline/hetAtr_indvs .
gzip *.vcf

# get a VCF of all species individually - tbh only need hetAtr separate at this point bc the outgroups are only single indv until we can get the 10x data to work with GATK
bcftools view -O z -s stiNae_male Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz
bcftools view -O z -s netAur_male Combined_hardFiltered.vcf > netAur.vcf.gz
bcftools view -O z -s oxyJam_male Combined_hardFiltered.vcf > oxyJam.vcf.gz

# output relevant stats - relatedness (unadjusted Ajk stat), heterozygosity on a per-individual basis, allele frequency at each site, nucleotide divergency on a per-site basis in 100kb windows, and a 012 genotype matrix 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.statsRel --relatedness 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.statsHet --het
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.statsFreq --freq 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.statsPi --window-pi 100000 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats012 --012 
vcftools --vcf Combined_hardFiltered.vcf --out allDucks.stats --012

mkdir qc
mv *stats* qc/

zgrep -v '\*' hetAtr.vcf.gz > hetAtr.clean.vcf
plink --vcf hetAtr.clean.vcf --make-bed --out hetAtr --allow-extra-chr
plink --bfile hetAtr --indep-pairwise 500 10 0.1 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --make-bed --extract hetAtr.prune.in --out hetAtr.ld_pruned --allow-extra-chr
plink --bfile hetAtr.ld_pruned --ibc --out hetAtr --allow-extra-chr

echo "x <- read.table(\"hetAtr.ibc\",header=T)" > hetAtr.plot.r
echo "pdf(\"hetAtr.ibc.pdf\",height=5,width=5)" >> hetAtr.plot.r
echo "plot(x\$Fhat1,x\$Fhat2,xlab=\"Fhat1\",ylab=\"Fhat2\")" >> hetAtr.plot.r
echo "dev.off()" >> hetAtr.plot.r
Rscript hetAtr.plot.r

plink --bfile hetAtr.ld_pruned --pca --out hetAtr --allow-extra-chr




