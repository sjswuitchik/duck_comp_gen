# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/vcfs/

module load bcftools/1.5-fasrc02 vcftools/0.1.14-fasrc01

cp ../shortRead_mapping_variantCalling/gatk/Combined_hardFiltered.vcf .
cp ~/github_repos/duck_comp_gen/NCBI_run/05_CompPopGen/02_MKpipeline/hetAtr_indvs .
gzip *.vcf

# get a VCF of all species individually - tbh only need hetAtr separate at this point bc the outgroups are only single indv until we can get the 10x data to work with GATK
bcftools view -O z -s stiNae_male Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz
bcftools view -O z -s netAur_male Combined_hardFiltered.vcf > netAur.vcf.gz
bcftools view -O z -s oxyJam_male Combined_hardFiltered.vcf > oxyJam.vcf.gz

# output relevant stats - relatedness (unadjusted Ajk stat), heterozygosity on a per-individual basis, allele frequency at each site, nucleotide divergency on a per-site basis in 100kb windows, and a 012 genotype matrix 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats --relatedness 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats --het
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats --freq 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats --window-pi 100000 
vcftools --gzvcf hetAtr.vcf.gz --out hetAtr.stats --012 
