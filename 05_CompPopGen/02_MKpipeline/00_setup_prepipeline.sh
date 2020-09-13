# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline, copy over these files: 
# hetAtr.vcf.gz
# stiNae.vcf.gz
# hetAtr_all_all_missingness_info.txt
# stiNae_all_all_missingness_info.txt

#### still need to figure out what the coverage output from snakemake will look like

cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_hetAtr.sorted.gtf .

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt
cut -f5 assembly_report.txt > scaffolds_only.txt
mv snpEff/data/hetAtr/sequences.fa .
grep -v -f scaffolds_only.txt sequences.fa > chr.only.seq.fa
### pick up here on Monday 

awk 'NR > 33 { print }' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | awk '{print $1, $5}' - > acckey
./replace_chrs.pl acckey cleaned_reordered_hetAtr.sorted.gtf > hetAtr.repchr.gtf

# download SnpEff
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
rm -r clinEff/
cd snpEff/

export R_LIBS_USER=$HOME/swuitchik/apps/R_3.6.1
export PATHW=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline
export PATHS=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline/snpEff
export INSHORT=hetAtr
export OUTGROUP=stiNae
# need to figure out if there should be an INLONG and OUTLONG...

mkdir -p data/$INSHORT
cd data/$INSHORT
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
gunzip GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz


mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna sequences.fa
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline/hetAtr.repchr.gtf .
mv hetAtr.repchr.gtf genes.gtf

cd ../..
 
conda activate mk_v2

java -jar $PATHS/snpEff.jar build -gtf22 -v $INSHORT
