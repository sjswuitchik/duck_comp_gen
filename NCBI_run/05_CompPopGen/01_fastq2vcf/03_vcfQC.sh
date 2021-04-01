# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/reseq_vcfs/

#conda create -n vcfqc -c bioconda plink vcftools bcftools r-base r-tidyverse admixture perl
source activate vcfqc

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
cd qc

# contstruct the 012 matrix with some wrangling for PCA, adegenet
cut -f2- hetAtr.stats012.012 > hetAtr.clean.012
paste hetAtr.stats012.012.indv hetAtr.clean.012 > hetAtr.matrix

# plink
zgrep -v '\*' hetAtr.vcf.gz > hetAtr.clean.vcf
plink --vcf hetAtr.clean.vcf --make-bed --out hetAtr --allow-extra-chr
plink --bfile hetAtr --indep-pairwise 500 50 0.1 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --make-bed --extract hetAtr.prune.in --out hetAtr.ld_pruned --allow-extra-chr --geno 0.95
plink --bfile hetAtr.ld_pruned --ibc --out hetAtr --allow-extra-chr

# plot of inbreeding coeficicients 
echo "x <- read.table(\"hetAtr.ibc\",header=T)" > hetAtr.plot.r
echo "pdf(\"hetAtr.ibc.pdf\",height=5,width=5)" >> hetAtr.plot.r
echo "plot(x\$Fhat1,x\$Fhat2,xlab=\"Fhat1\",ylab=\"Fhat2\")" >> hetAtr.plot.r
echo "dev.off()" >> hetAtr.plot.r
Rscript hetAtr.plot.r

# PCA 
plink --bfile hetAtr.ld_pruned --pca --out hetAtr --allow-extra-chr

echo "pdf(\"hetAtr.pca.pdf\",height=8,width=5)" > hetAtr.pca.plot
echo "par(mfrow=c(2,1),mar=c(4,4,2,2))" >> hetAtr.pca.plot
echo "d <- read.table(\"hetAtr.eigenval\")" >> hetAtr.pca.plot
echo "plot(c(seq(1,length(d\$V1),by=1)),d\$V1/sum(d\$V1)*100,xlab=\"PC\",ylab=\"Percent Variance Explained\")" >> hetAtr.pca.plot
echo "d <- read.table(\"hetAtr.eigenvec\")" >> hetAtr.pca.plot
echo "plot(d\$V3,d\$V4,cex=0.5,xlab=\"PC 1\",ylab = \"PC 2\")" >> hetAtr.pca.plot
echo "dev.off()" >> hetAtr.pca.plot
Rscript hetAtr.pca.plot

# admixture

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt
sed 's/\r$//g' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | grep -v "^#" | cut -f1,3,5 > hetAtr_chr_key
awk '{print $3, $2}' hetAtr_chr_key > acckey
# manually editing the chrom names in acckey to include 'chr' for now, see if sed or something similar will work later; chr 34 = W, chr 35 = Z, chr 36 = MT
# nb: manually editing that file sucked, figure out a better way to do this
sed -i 's/na/\d0/g' acckey 
./replace_chrs.pl acckey hetAtr.ld_pruned.bim > hetAtr.repl.ld_pruned.bim
mv hetAtr.repl.ld_pruned.bim hetAtr.ld_pruned.bim

for K in {2..5}
do
	admixture --cv hetAtr.ld_pruned.bed $K > hetAtr.${K}.admix.log 2> hetAtr.${K}.admix.err
done

cat hetAtr.*.admix.log | grep CV | perl -pi -e 's/.+=//' | perl -pi -e 's/\): /\t/' > hetAtr.CV

# plot admixture output
echo "x <- read.table(\"hetAtr.CV\")" > hetAtr.admixture.plot.r
echo "pdf(\"hetAtr.admix.pdf\",height=10,width=5)" >> hetAtr.admixture.plot.r
echo "par(mfrow=c(6,1),mar=c(4,4,2,2))" >> hetAtr.admixture.plot.r
echo "plot(x\$V1,x\$V2,xlab=\"K\",ylab=\"CV\")" >> hetAtr.admixture.plot.r
echo "Q1 <- as.matrix(read.table(\"hetAtr.ld_pruned.1.Q\"))" >> hetAtr.admixture.plot.r
echo "Q2 <- as.matrix(read.table(\"hetAtr.ld_pruned.2.Q\"))" >> hetAtr.admixture.plot.r
echo "Q3 <- as.matrix(read.table(\"hetAtr.ld_pruned.3.Q\"))" >> hetAtr.admixture.plot.r
echo "Q4 <- as.matrix(read.table(\"hetAtr.ld_pruned.4.Q\"))" >> hetAtr.admixture.plot.r
echo "Q5 <- as.matrix(read.table(\"hetAtr.ld_pruned.5.Q\"))" >> hetAtr.admixture.plot.r
echo "barplot(t(Q1),col=rainbow(6))" >> hetAtr.admixture.plot.r
echo "barplot(t(Q2),col=rainbow(2))" >> hetAtr.admixture.plot.r
echo "barplot(t(Q3),col=rainbow(3))" >>	hetAtr.admixture.plot.r
echo "barplot(t(Q4),col=rainbow(4))" >>	hetAtr.admixture.plot.r
echo "barplot(t(Q5),col=rainbow(5))" >> hetAtr.admixture.plot.r
echo "dev.off()" >> hetAtr.admixture.plot.r

Rscript hetAtr.admixture.plot.r




## VariantQC
module load jdk/10.0.1-fasrc01 htslib/1.5-fasrc02 bcftools/1.5-fasrc02
cd /n/holyscratch01/informatics/swuitchik/ducks/snakemake/reseq_vcfs
cp hetAtr.vcf.gz stiNae.vcf.gz ../hetAtr_run/genome/hetAtr.*  ../hetAtr_stiNae_qc
cd ../hetAtr_stiNae_qc

# sort
source activate gatk
picard SortVcf -Xmx8g -I hetAtr.vcf.gz -O hetAtr.sorted.vcf.gz
picard SortVcf -Xmx8g -I stiNae.vcf.gz -O stiNae.sorted.vcf.gz

# index
gatk IndexFeatureFile -I hetAtr.sorted.vcf.gz
gatk IndexFeatureFile -I stiNae.sorted.vcf.gz

# get VariantQC tool & run
wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.24/DISCVRSeq-1.24.jar
mv DISCVRSeq-1.24.jar DISCVRSeq.jar
java -jar DISCVRSeq.jar VariantQC -R hetAtr.fa -V hetAtr.sorted.vcf.gz -O hetAtr.html
java -jar DISCVRSeq.jar VariantQC -R hetAtr.fa -V stiNae.sorted.vcf.gz -O stiNae.html

## try VariantQC with updated Filter VCF 
cp ../hetAtr_run/Combined_hardFiltered_updatedFilter.vcf.gz* ../hetAtr_run/hetAtr_indvs .
bcftools view -O z -S hetAtr_indvs -G -a Combined_hardFiltered_updatedFilter.vcf.gz > hetAtr.filtered.vcf.gz
bcftools view -O z -s stiNae_male -G -a Combined_hardFiltered_updatedFilter.vcf.gz > stiNae.filtered.vcf.gz

# sort
picard SortVcf -Xmx8g -I hetAtr.filtered.vcf.gz -O hetAtr.filtered.sorted.vcf.gz
picard SortVcf -Xmx8g -I stiNae.filtered.vcf.gz -O stiNae.filtered.sorted.vcf.gz

# index
gatk IndexFeatureFile -I hetAtr.filtered.sorted.vcf.gz
gatk IndexFeatureFile -I stiNae.filtered.sorted.vcf.gz

# run VariantQC
java -jar DISCVRSeq.jar VariantQC -R hetAtr.fa -V hetAtr.filtered.sorted.vcf.gz -O hetAtr.filtered.html
java -jar DISCVRSeq.jar VariantQC -R hetAtr.fa -V stiNae.filtered.sorted.vcf.gz -O stiNae.filtered.html

# output stats
vcftools --gzvcf hetAtr.filtered.sorted.vcf.gz --out hetAtr.rel --relatedness2
vcftools --gzvcf hetAtr.filtered.sorted.vcf.gz --out hetAtr.10kb --TajimaD 10000
vcftools --gzvcf hetAtr.filtered.sorted.vcf.gz --out hetAtr.statsPi --window-pi 100000 

zgrep -v '\*' hetAtr.filtered.sorted.vcf.gz > hetAtr.filtered.sorted.clean.vcf
plink --vcf hetAtr.filtered.sorted.clean.vcf --make-bed --out hetAtr --allow-extra-chr
plink --bfile hetAtr --indep-pairwise 500 50 0.1 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --make-bed --extract hetAtr.prune.in --out hetAtr.ld_pruned --allow-extra-chr
plink --bfile hetAtr.ld_pruned --ibc --out hetAtr --allow-extra-chr







