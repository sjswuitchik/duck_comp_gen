# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/gatk/

# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline

module load bcftools/1.5-fasrc02 bedtools2/2.26.0-fasrc01 perl/5.26.1-fasrc01 python/2.7.14-fasrc02

cp Combined_hardFiltered.vcf missing_data_per_ind.txt /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline 
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline

# separate hetAtr and stiNae VCFs
bcftools view -O z -s stiNae_ind01 Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz

# separate hetAtr and stiNae missingness data
grep -v stiNae_ind01 missing_data_per_ind.txt > hetAtr_missing_data.txt
sed -n 1p missing_data_per_ind.txt > stiNae_missing_data.txt
grep -v -f hetAtr_missing_data.txt missing_data_per_ind.txt > temp
cat temp >> stiNae_missing_data.txt
rm temp

# calculate coverage 
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/fastq2bam/01_mappedReads

for file in *_dedup.bam;
do
  bedtools genomecov -bga -ibam $file -g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/hetAtr/genome/hetAtr.fa > $file.statscov.bg
done

wget https://github.com/shenwei356/brename/releases/download/v2.10.0/brename_linux_amd64.tar.gz
tar zxvf brename_linux_amd64.tar.gz 
rm brename_linux_amd64.tar.gz 
chmod +x ./brename

./brename -p "_dedup.bam." -r "." -R 

cp *.statscov.bg /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/coverage_stats
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/coverage_stats

gzip stiNae_ind01.statscov.bg

sbatch unioncov.sh

gzip -dc hetAtr_union.bg.gz | ./sum_cov.awk
gunzip coverage*.gz
sed '1d' coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_clean_merged.bed
sed '1d' coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_low_merged.bed
sed '1d' coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_high_merged.bed

gzip -dc stiNae_ind01.statscov.bg.gz | ./sum_cov.awk
gunzip coverage*.gz
sed '1d' coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_clean_merged.bed
sed '1d' coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_low_merged.bed
sed '1d' coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_high_merged.bed

cp -v hetAtr_coverage_sites_clean_merged.bed stiNae_coverage_sites_clean_merged.bed ..
cd ..

# using GTF from CESAR annotation for snpEff database build
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/snpEff/data/hetAtr
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03b_cesar/output_gtfs/cleaned_reordered_hetAtr.sorted.gtf .
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt
sed 's/\r$//g' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | grep -v "^#" | cut -f1,5 > hetAtr_key
./replace_chrs.pl hetAtr_key cleaned_reordered_hetAtr.sorted.gtf > genes.gtf
python2 gtf2gff.py genes.gtf > genes.gff
gzip genes.gff

cd ../..
module load Anaconda/5.0.1-fasrc02
source activate mk_v2
java -jar snpEff.jar build -gff3 -v hetAtr









