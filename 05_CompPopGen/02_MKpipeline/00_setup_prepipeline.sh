# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline

module load Anaconda3/2019.10
#conda create -c bioconda -c conda-forge -n busco busco=4.0.6 gffread

# only using previously build busco env for the gffread util
conda activate busco

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline, copy over these files: 
# hetAtr.vcf.gz
# stiNae.vcf.gz
# hetAtr_all_all_missingness_info.txt
# stiNae_all_all_missingness_info.txt

#### still need to figure out what the coverage output from snakemake will look like

mkdir gffs
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_hetAtr.sorted.gtf gffs/
cd gffs/

gffread cleaned_reordered_hetAtr.sorted.gtf -o hetAtr.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt

awk 'NR > 33 { print }' GCA_011075105.1_BPBGC_Hatr_1.0_assembly_report.txt | awk '{print $1, $5}' - > acckey
./replace_chrs.pl acckey hetAtr.gff > hetAtr.repchr.gff

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
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline/gffs/hetAtr.repchr.gff .
mv hetAtr.gff genes.gff
gzip genes.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz sequences.fa.gz

cd ../..

conda deactivate 
conda activate mk_v2

java -jar $PATHS/snpEff.jar build -gff3 -v $INSHORT
