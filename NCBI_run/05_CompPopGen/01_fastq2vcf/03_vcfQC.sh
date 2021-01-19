# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/vcfs/

module load bcftools/1.5-fasrc02

cp ../shortRead_mapping_variantCalling/gatk/Combined_hardFiltered.vcf .
cp ~/github_repos/duck_comp_gen/NCBI_run/05_CompPopGen/02_MKpipeline/hetAtr_indvs .

bcftools view -O z -s stiNae_male Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz
bcftools view -O z -s netAur_male Combined_hardFiltered.vcf > netAur.vcf.gz
bcftools view -O z -s oxyJam_male Combined_hardFiltered.vcf > oxyJam.vcf.gz
