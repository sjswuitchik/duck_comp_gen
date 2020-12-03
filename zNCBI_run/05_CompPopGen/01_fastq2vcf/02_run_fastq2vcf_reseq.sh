# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling

module load Anaconda/5.0.1-fasrc01 bcftools/1.5-fasrc02
#conda create --name snakemake -c bioconda snakemake samtools picard

# make directories for data
for file in netAur oxyJam stiNae;
do
  mkdir -p data/$file/fastqs
  mkdir -p data/$file/genome
done

# make a v2 directory for hetAtr to include both the reseq data and the genome data (was not included in initial run with REU data)
mkdir -p data/hetAtr_v2/fastqs
mkdir -p data/hetAtr_v2/genome

# copy fastqs over 
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/BGI_duck_reseq/CNB0005/* /n/holylfs/LABS/informatics/swuitchik/ducks/raw_fastqs/01_stiNae/* data/stiNae/fastqs
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/BGI_duck_reseq/CNB0006/* /n/holylfs/LABS/informatics/swuitchik/ducks/raw_fastqs/03_netAur/* data/netAur/fastqs
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/BGI_duck_reseq/CNB0007/* /n/holylfs/LABS/informatics/swuitchik/ducks/raw_fastqs/02_oxyJam/* data/oxyJam/fastqs
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/hetAtr/fastqs/hetAtr* /n/holylfs/LABS/informatics/swuitchik/ducks/raw_fastqs/04_hetAtr/* data/hetAtr_v2/fastqs

# get reference genomes
cd data/hetAtr_v2/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz hetAtr.fa.gz

cd ../../oxyJam/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/077/185/GCF_011077185.1_BPBGC_Ojam_1.0/GCF_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
mv GCF_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.fa.gz

cd ../../netAur/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz netAur.fa.gz

cd ../../stiNae/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz stiNae.fa.gz

cd ../..
for file in netAur oxyJam stiNae hetAtr;
do
  gunzip $file/genome/*.gz
done

# concat fastqs (in same order!)
./concat_fastqs_reseq.sh 

# move original fastqs so only the prepped fastqs are called in the pipeline
for file in hetAtr_v2 oxyJam netAur stiNae;
do
  cd $file/fastqs
  mkdir orig_fastqs
done

cd stiNae/fastqs
mv stiNae_female_L00* orig_fastqs/
cd ../../oxyJam/fastqs
mv oxyJam_female_L00* orig_fastqs/
cd ../../netAur/fastqs
mv netAur_female_L00* orig_fastqs/
cd ../../hetAtr_v2/fastqs
mv DGAB-CNB0004-CN4-lib1_S* orig_fastqs/
cd ../../

#### while working out the issues with 10x data memory & Picard MarkDuplicates, aligning all non-10x data to black-headed duck genome ####
mkdir -p allDucks/fastqs
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/hetAtr/fastqs/hetAtr* allDucks/fastqs/
for file in oxyJam netAur stiNae;
do
  cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/01_fastq2vcf/shortRead_mapping_variantCalling/data/$file/fastqs/$file_male* allDucks/fastqs/
done

cd ../..

# move all previous run data for hetAtr& stiNae MK pipeline input so it doesn't get overwritten
mkdir hetAtr_stiNae_MK
cd shortRead_mapping_variantCalling/
mv fastq2bam/ intervalFiles/ gatk/ ../hetAtr_stiNae_MK

# run allDucks
sbatch run_fastq2bam.sh
sbatch run_intervals.sh
sbatch run_bam2vcf_gatk.sh

