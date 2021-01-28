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

# plink 
zgrep -v '\*' hetAtr.vcf.gz > hetAtr.clean.vcf
plink --vcf hetAtr.clean.vcf --make-bed --out hetAtr --allow-extra-chr
plink --bfile hetAtr --indep-pairwise 500 10 0.1 --out hetAtr --allow-extra-chr
plink --bfile hetAtr --make-bed --extract hetAtr.prune.in --out hetAtr.ld_pruned --allow-extra-chr
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
for K in {2..5}
do
	admixture --cv ${SPECIES}.ld_pruned.bed $K > ${SPECIES}.${K}.admix.log 2> ${SPECIES}.${K}.admix.err
done

cat hetAtr.*.admix.log | grep CV | perl -pi -e 's/.+=//' | perl -pi -e 's/\): /\t/' > hetAtr.CV

# plot admixture output
echo "x <- read.table(\"hetAtr.CV\")" > hetAtr.admixture.plot.r
echo "pdf(\"hetAtr.admix.pdf\",height=10,width=5)" >> hetAtr.admixture.plot.r
echo "par(mfrow=c(5,1),mar=c(4,4,2,2))" >> hetAtr.admixture.plot.r
echo "plot(x\$V1,x\$V2,xlab=\"K\",ylab=\"CV\")" >> hetAtr.admixture.plot.r
echo "Q2 <- as.matrix(read.table(\"hetAtr.ld_pruned.2.Q\"))" >> hetAtr.admixture.plot.r
echo "Q3 <- as.matrix(read.table(\"hetAtr.ld_pruned.3.Q\"))" >> hetAtr.admixture.plot.r
echo "Q4 <- as.matrix(read.table(\"hetAtr.ld_pruned.4.Q\"))" >> hetAtr.admixture.plot.r
echo "Q5 <- as.matrix(read.table(\"hetAtr.ld_pruned.5.Q\"))" >> hetAtr.admixture.plot.r
echo "barplot(t(Q2),col=rainbow(2))" >> hetAtr.admixture.plot.r
echo "barplot(t(Q3),col=rainbow(3))" >>	hetAtr.admixture.plot.r
echo "barplot(t(Q4),col=rainbow(4))" >>	hetAtr.admixture.plot.r
echo "barplot(t(Q5),col=rainbow(5))" >> hetAtr.admixture.plot.r
echo "dev.off()" >> hetAtr.admixture.plot.r

Rscript hetAtr.admixture.plot.r

