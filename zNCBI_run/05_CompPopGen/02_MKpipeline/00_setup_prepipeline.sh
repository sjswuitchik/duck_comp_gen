# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/gatk/

# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline

module load bcftools/1.5-fasrc02

cp Combined_hardFiltered.vcf missing_data_per_ind.txt /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline 
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline

# pull out the stiNae individual to its own VCF
bcftools view -O z -s stiNae_ind01 Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz

# pull out stiNae from missingness file & rename
grep -v stiNae_ind01 missing_data_per_ind.txt > hetAtr_missing_data.txt
sed -n 1p missing_data_per_ind.txt > stiNae_missing_data.txt
grep -v -f hetAtr_missing_data.txt missing_data_per_ind.txt > temp
cat temp >> stiNae_missing_data.txt
rm temp

