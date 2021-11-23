# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/

module load Anaconda/5.0.1-fasrc01 bcftools/1.5-fasrc02
#conda create --name snakemake -c bioconda snakemake samtools picard

# clone pipeline & repo 
git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling
cd shortRead_mapping_variantCalling

# make directories for data
mkdir -p data/hetAtr/fastqs
mkdir -p data/hetAtr/genome

# copy fastqs over 
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen_old/orig_fastqs/stiNae_male/*.gz /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/hetAtr_reseq/*.gz data/hetAtr/fastqs

# get reference genome
cd data/hetAtr/genome
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/hetAtr/hetAtr.fasta .
mv hetAtr.fasta hetAtr.fa

# concat fastqs (in same order!) for individuals 
cd data/hetAtr/fastqs
./concat_fastqs.sh 
mkdir orig_fastqs/
mv hetAtr*_L*.fastq.gz orig_fastqs/

sbatch run_fastq2bam.sh
sbatch run_intervals.sh
sbatch run_bam2vcf_gatk.sh
