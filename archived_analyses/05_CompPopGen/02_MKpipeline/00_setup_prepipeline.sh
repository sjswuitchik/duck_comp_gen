# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline

module load Anaconda3/2019.10 samtools/1.10-fasrc01

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline, copy over these files: 
# hetAtr.vcf.gz
# stiNae.vcf.gz
# hetAtr_missing_data.txt
# stiNae_missing_data.txt
# hetAtr_coverage_sites_clean_merged.bed
# stiNae_coverage_sites_clean_merged.bed

# download SnpEff
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
rm -r clinEff/
cd snpEff/
mkdir -p data/$INSHORT
cd data/$INSHORT
cd ../..

export R_LIBS_USER=$HOME/swuitchik/apps/R_3.6.1
export PATHW=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline
export PATHS=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline/snpEff
export INSHORT=hetAtr
export OUTSHORT=stiNae


# download genome and pull out only chromosomes (to match chromosome only WGA annotation)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
gunzip GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
samtools faidx -r chr_only.txt GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna > chr.only.seq.fa
mv chr.only.seq.fa sequences.fa
rm GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.fai
mv sequences.fa snpEff/data/hetAtr/

# get assembly report to replace chromosome names in GTF (from het1 -> CM021731.1, etc.) 
conda activate mk_v2
conda install -c bioconda agat # needed for gtf conversion

cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_hetAtr.sorted.gtf .
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt
awk 'NR > 33 { print }' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | awk '{print $1, $5}' - > acckey
./replace_chrs.pl acckey cleaned_reordered_hetAtr.sorted.gtf > hetAtr.repchr.gtf

agat_convert_sp_gxf2gxf.pl -g hetAtr.repchr.gtf --gvo 3 -o hetAtr.repchr.gff

mv hetAtr.repchr.gff genes.gff
mv genes.gff snpEff/data/hetAtr

# build snpEff database
java -jar $PATHS/snpEff.jar build -gff3 -v $INSHORT
