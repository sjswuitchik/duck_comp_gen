# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run
# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline. Copy or move over the VCF and the missing data files from the snakemake pipeline output and copy genome dir from fastq2bam as well 

module load bcftools/1.5-fasrc02 bedtools2/2.26.0-fasrc01 perl/5.26.1-fasrc01 Anaconda3/2020.11

# update filtering of VCF and clean up header
sbatch run_gatkUpdate.sh Combined_hardFiltered hetAtr
./vcfHeaderCleanup 

# separate hetAtr and stiNae VCFs
bcftools view -O z -s stiNae_male Combined_hardFiltered_updatedFilter.vcf.gz > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered_updatedFilter.vcf.gz > hetAtr.vcf.gz

# separate hetAtr and stiNae missingness data
sed -n 1p missing_data_per_ind.txt > hetAtr_missing_data.txt
grep -f hetAtr_indvs missing_data_per_ind.txt > temp
cat temp >> hetAtr_missing_data.txt

sed -n 1p missing_data_per_ind.txt > stiNae_missing_data.txt
grep stiNae_male missing_data_per_ind.txt > temp
cat temp >> stiNae_missing_data.txt
rm temp

# calculate coverage 
# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/fastq2bam_hetAtr/01_mappedReads
for file in *_dedup.bam;
do
  bedtools genomecov -bga -ibam $file -g ../../genome/hetAtr.fa > $file.bg
done

mkdir /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/coverage
mv *.bg /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/coverage

cd ../..
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo

./faToTwoBit genome/hetAtr.fa hetAtr.2bit
./twoBitInfo hetAtr.2bit stdout | sort -k2rn > hetAtr.chrom.sizes

cd coverage/
# nb: to download binaries: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/
sbatch run_bedg2bw.sh 

ls hetAtr*.bw > hetAtr_list
/n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/./bigWigMerge -inList hetAtr_list hetAtr.merge.bg
/n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/./bedGraphToBigWig hetAtr.merge.bg ../hetAtr.chrom.sizes hetAtr.merge.bw
awk 'BEGIN{FS=OFS="\t"}{print $1, 0, $2, $1}' ../hetAtr.chrom.sizes > ../hetAtr.genome.bed
/n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/./bigWigAverageOverBed hetAtr.merge.bw ../hetAtr.genome.bed hetAtr.summary.tab 
/n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage/./bigWigAverageOverBed stiNae_male.bw ../hetAtr.genome.bed stiNae.summary.tab

sbatch write_coverage_beds.sh hetAtr
mean=$(awk '{sum = sum+$4}{size=size+$2}{avg=sum/size}END{print avg}' stiNae.summary.tab)
gzip -dc stiNae_male_dedup.bam.bg.sorted.gz | awk -v avg="$mean" -v spp=stiNae -f ../sum_cov.awk

sed '1d' hetAtr_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_clean_merged.bed
sed '1d' hetAtr_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_low_merged.bed
sed '1d' hetAtr_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_high_merged.bed

sed '1d' stiNae_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_clean_merged.bed
sed '1d' stiNae_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_low_merged.bed
sed '1d' stiNae_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_high_merged.bed

cp -v hetAtr_coverage_sites_clean_merged.bed stiNae_coverage_sites_clean_merged.bed /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline
cd /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline

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
#conda create -n r -c bioconda r-base r-tidyverse
source activate r
Rscript reformat.R
# add a gene ID to the GFF with chicken-based genes from translation files
./duck2chick-gtf.sh hetAtr.gff > hetAtr_final.gtf
# rename for snpEff
cp hetAtr_final.gtf genes.gtf
cp genes.gtf ..
cd ..
cp genes.gtf snpEff/data/hetAtr
gzip snpEff/data/hetAtr/genes.gtf
conda deactivate



source activate mk
cd snpEff
java -jar snpEff.jar build -gff3 -v hetAtr









