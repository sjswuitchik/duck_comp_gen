# in /n/holyscratch01/informatics/swuitchik/ducks/
# there is a little bit of prep that needs to be done before the output from the snakemake pipeline will be suitable to work in the MK pipeline. Copy or move over the VCF and the missing data files from the snakemake pipeline output and copy genome dir from fastq2bam as well 

module load bcftools/1.5-fasrc02

# clone MK repo 
git clone https://github.com/sjswuitchik/compPopGen_ms.git
mv compPopGen_ms/MKpipline .
rm -rf compPopGen_ms

# update filtering of VCF and clean up header
cd /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run
sbatch run_gatkUpdate.sh Combined_hardFiltered hetAtr
./vcfHeaderCleanup 

# separate hetAtr and stiNae VCFs
bcftools view -O z -s stiNae_male Combined_hardFiltered_updatedFilter.vcf.gz > stiNae.vcf.gz
bcftools view -O z -S hetAtr_indvs Combined_hardFiltered_updatedFilter.vcf.gz > hetAtr.vcf.gz

# copy VCFs to MK working dir
cp -v stiNae.vcf.gz hetAtr.vcf.gz /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline

# separate hetAtr and stiNae missingness data
sed -n 1p missing_data_per_ind.txt > hetAtr_missing_data.txt
grep -f hetAtr_indvs missing_data_per_ind.txt > temp
cat temp >> hetAtr_missing_data.txt

sed -n 1p missing_data_per_ind.txt > stiNae_missing_data.txt
grep stiNae_male missing_data_per_ind.txt > temp
cat temp >> stiNae_missing_data.txt
rm temp

# copy missing data files to MK working dir
cp -v stiNae_missing_data.txt hetAtr_missing_data.txt /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline

# calculate coverage 
# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/fastq2bam_hetAtr/01_mappedReads
#conda create -n coverage -c bioconda ucsc-fatotwobit ucsc-twobitinfo bedtools ucsc-bedgraphtobigwig ucsc-bigwigmerge ucsc-bigwigaverageoverbed 
source activate coverage
# create bedgraphs
for file in *_dedup.bam;
do
  bedtools genomecov -bga -ibam $file -g ../../genome/hetAtr.fa > $file.bg
done

mkdir /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/coverage
mv *.bg /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_run/coverage

# create chromosome sizes file
cd ../..
faToTwoBit genome/hetAtr.fa hetAtr.2bit
twoBitInfo hetAtr.2bit stdout | sort -k2rn > hetAtr.chrom.sizes

cd coverage/
# nb: to download binaries: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/
# convert bedgraphs to bigwigs 
sbatch run_bedg2bw.sh 

# merge bigwigs
ls hetAtr*.bw > hetAtr_list
bigWigMerge -inList hetAtr_list hetAtr.merge.bg

# convert merged bedgraph back to bigwig
bedGraphToBigWig hetAtr.merge.bg ../hetAtr.chrom.sizes hetAtr.merge.bw

# create chrom sizes bed file
awk 'BEGIN{FS=OFS="\t"}{print $1, 0, $2, $1}' ../hetAtr.chrom.sizes > ../hetAtr.genome.bed

# summarize merged bigwigs
bigWigAverageOverBed hetAtr.merge.bw ../hetAtr.genome.bed hetAtr.summary.tab 
bigWigAverageOverBed stiNae_male.bw ../hetAtr.genome.bed stiNae.summary.tab

# write out coverage bed files (low, high, clean sites)
gzip hetAtr.merge.bg
sbatch write_coverage_beds.sh hetAtr

mean=$(awk '{sum = sum+$4}{size=size+$2}{avg=sum/size}END{print avg}' stiNae.summary.tab)
gzip -dc stiNae_male_dedup.bam.bg.sorted.gz | awk -v avg="$mean" -v spp=stiNae -f ../sum_cov.awk

# sort & merge coverage bed files
sed '1d' hetAtr_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_clean_merged.bed
sed '1d' hetAtr_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_low_merged.bed
sed '1d' hetAtr_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > hetAtr_coverage_sites_high_merged.bed

sed '1d' stiNae_coverage_sites_clean.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_clean_merged.bed
sed '1d' stiNae_coverage_sites_low.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_low_merged.bed
sed '1d' stiNae_coverage_sites_high.bed | bedtools sort -i - | bedtools merge -i - > stiNae_coverage_sites_high_merged.bed

conda deactivate 

# copy files to working dir for MK pipeline
cp -v hetAtr_coverage_sites_clean_merged.bed stiNae_coverage_sites_clean_merged.bed /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline
cd /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline

# using GTF from Comparative Augustus annotation for snpEff database build - CompAug GTF needs to be translated using OrthoFinder orthologue comparison of hetAtr v galGal
# in /n/holyscratch01/informatics/swuitchik/ducks/MKpipeline/GFFtranslate
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/Orthologues/Orthologues_hetAtr.translated/hetAtr.translated__v__galGal.translated.tsv .
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/hetAtr.gff .

# get galGal6 genome annotation
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

# add a gene ID to the GTF with chicken-based genes from translation files
./duck2chick-gtf.sh hetAtr.gff > hetAtr_final.gtf

# rename for snpEff
cp hetAtr_final.gtf genes.gtf
cp genes.gtf ..
cd ..
cp genes.gtf snpEff/data/hetAtr
gzip snpEff/data/hetAtr/genes.gtf
conda deactivate

# grab hetAtr genome for snpEff
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/genomes/hetAtr.ncbi.fasta snpEff/data/hetAtr/
mv snpEff/data/hetAtr/hetAtr.ncbi.fasta sequences.fa
gzip snpEff/data/hetAtr/sequences.fa

# build snpEff database 
source activate mk
cd snpEff
snpEff -Xmx8g build -c snpEff.config -gtf22 -v hetAtr
