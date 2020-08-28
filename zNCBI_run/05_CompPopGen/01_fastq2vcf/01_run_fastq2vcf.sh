# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/

module load Anaconda/5.0.1-fasrc01 bcftools/1.5-fasrc02
#conda create --name snakemake -c bioconda snakemake 

# clone pipeline & repo 
git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling
cd shortRead_mapping_variantCalling

# make directories for data
mkdir -p data/hetAtr/fastqs
mkdir -p data/hetAtr/genome

# copy fastqs over 
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/stiNae_male/*.gz /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/hetAtr_reseq/*.gz data/hetAtr/fastqs

# get reference genome
cd data/hetAtr/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
gunzip GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna hetAtr.fa

# concat fastqs (in same order!) for individuals 
cd data/hetAtr/fastqs
./concat_fastqs.sh 
mkdir orig_fastqs/
mv hetAtr*_L*.fastq.gz orig_fastqs/

sbatch run_fastq2bam.sh
sbatch run_bam2vcf_gatk.sh

# pull out the stiNae individual to its own VCF
bcftools view -c1 -Oz -s stiNae_ind01 -o stiNae.vcf.gz Combined_hardFiltered.vcf

mv Combined_hardFiltered.vcf hetAtr.vcf.gz


