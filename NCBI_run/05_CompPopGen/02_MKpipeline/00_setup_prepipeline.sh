# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling

# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline

module load bcftools/1.5-fasrc02 bedtools2/2.26.0-fasrc01 perl/5.26.1-fasrc01 Anaconda/5.0.1-fasrc02

cp gatk/Combined_hardFiltered.vcf manual/missing_data_per_ind.txt ../../02_MK_pipeline 
cd ../../02_MK_pipeline

# separate hetAtr and stiNae VCFs
bcftools view -O z -s stiNae_male Combined_hardFiltered.vcf > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered.vcf > hetAtr.vcf.gz

# separate hetAtr and stiNae missingness data
sed -n 1p missing_data_per_ind.txt > hetAtr_missing_data.txt
grep -f hetAtr_indvs missing_data_per_ind.txt > temp
cat temp >> hetAtr_missing_data.txt

sed -n 1p missing_data_per_ind.txt > stiNae_missing_data.txt
grep stiNae_male missing_data_per_ind.txt > temp
cat temp >> stiNae_missing_data.txt
rm temp

# calculate coverage 
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/fastq2bam/01_mappedReads/
for file in *_dedup.bam;
do
  bedtools genomecov -bga -ibam $file -g ../../data/allDucks/genome/hetAtr.fa > $file.statscov.bg
done

wget https://github.com/shenwei356/brename/releases/download/v2.10.0/brename_linux_amd64.tar.gz
tar zxvf brename_linux_amd64.tar.gz 
rm brename_linux_amd64.tar.gz 
chmod +x ./brename

./brename -p "_dedup.bam." -r "." -R 

cp *.statscov.bg /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/coverage_stats
cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/coverage_stats

gzip stiNae_male.statscov.bg

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/allDucks/genome/hetAtr.fa .

./faToTwoBit hetAtr.fa hetAtr.2bit
./twoBitInfo hetAtr.2bit stdout | sort -k2rn > hetAtr.chrom.sizes

sbatch unioncov.sh

gzip -dc hetAtr_union.bg.gz | ./sum_cov.awk
sed '1d' hetAtr_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_clean_merged.bed
sed '1d' hetAtr_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_low_merged.bed
sed '1d' hetAtr_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_high_merged.bed

gzip -dc stiNae_male.statscov.bg.gz | ./stiNae_sum_cov.awk
sed '1d' stiNae_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_clean_merged.bed
sed '1d' stiNae_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_low_merged.bed
sed '1d' stiNae_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_high_merged.bed

cp -v hetAtr_coverage_sites_clean_merged.bed stiNae_coverage_sites_clean_merged.bed ..
cd ..

# using GFF from Comparative Augustus annotation for snpEff database build
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline/snpEff/data/hetAtr
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/hetAtr.ncbi.fasta .
mv hetAtr.ncbi.fasta sequences.fa
gzip sequences.fa
# CompAug GFF needs to be translated using OrthoFinder orthologue comparison of hetAtr v galGal
cd ../../../
mkdir hetAtr_translation
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/OrthoFinder_jan2021/run_ortho/OrthoFinder/Results_Feb01/Orthologues/Orthologues_hetAtr.translated/hetAtr.translated__v__galGal.translated.tsv hetAtr_translation
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/hetAtr.gff hetAtr_translation/
cd hetAtr_translation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz
# create translation file for galGal transcript to gene
grep -v '#' GCF_000002315.6_GRCg6a_genomic.gff | awk '{if ($3 == "mRNA") print $0;}' | python3 galGenes_trans.py > transGene.txt
# create translation file for hetAtr to galGal transcripts
sed '1d' hetAtr.translated__v__galGal.translated.tsv | cut -f2,3 > hetGal_trans.tsv
# quickly reformat transGene file
# conda create -n agat -c conda-forge -c bioconda agat r-base r-tidyerse 
source activate agat
Rscript reformat.R
# replace using AGAT scripts 
agat_sp_manage_IDs.pl 



mv hetAtr.gff genes.gff
gzip genes.gff


cd ../../..
module load Anaconda/5.0.1-fasrc02
source activate mk
java -jar snpEff.jar build -gff3 -v hetAtr









