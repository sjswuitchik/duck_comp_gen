# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling

# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline

# pull out the stiNae individual to its own VCF
bcftools view -c1 -Oz -s stiNae_ind01 -o stiNae.vcf.gz Combined_hardFiltered.vcf

mv Combined_hardFiltered.vcf hetAtr.vcf.gz

# move VCFs and missingness over to MK dir
mv stiNae.vcf.gz hetAtr.vcf.gz /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline
mv missing_data_per_ind.txt /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline

cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline

# pull out stiNae from missingness file & rename
grep -v stiNae_ind01 missing_data_per_ind.txt > hetAtr_all_all_missingness_info.txt
grep -v -f hetAtr_all_all_missingness_info.txt missing_data_per_ind.txt > stiNae_all_all_missingness_info.txt

