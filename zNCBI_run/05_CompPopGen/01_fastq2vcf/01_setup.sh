# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/

module load Anaconda/5.0.1-fasrc01
#conda create --name snakemake -c bioconda snakemake 

# clone pipeline & repo 
git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling
cd shortRead_mapping_variantCalling

# make directories for data
mkdir -p data/stiNae/fastqs
mkdir -p data/hetAtr/fastqs
mkdir -p data/hetAtr/genome

# copy fastqs over 
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/stiNae_male/*.gz data/stiNae/fastqs
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/hetAtr_reseq/*.gz data/hetAtr/fastqs

# get reference genome
cd data/hetAtr/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
gunzip GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna hetAtr.fa

# concat fastqs (in same order!) for individuals 
./data/hetAtr/fastqs/concat_fastqs.sh 

sbatch run_fastq2bam.sh


