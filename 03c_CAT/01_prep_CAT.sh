######################################
## Gene annotation with CAT - setup ##
######################################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cat
# also testing in /scratch/swuitchik on bioinf01 for optimization and debugging
# https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit

module load Anaconda3/2019.10 
#conda create -n gmap -c bioconda gmap samtools 

# extract only focal and reference species from WGA
mkdir input_data
cd input_data
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
hal2maf galloanserae.hal galloForCAT.maf --refGenome galGal --noAncestors --noDupes --targetGenomes galGal,hetAtr,netAur,oxyJam,stiNae
maf2hal galloForCAT.maf galloForCAT.hal --refGenome galGal
exit 

# get filtered GFF for galGal
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.filtpy.gff . 

# get galGal protein alignment
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_protein.faa.gz

# get indexed fastas for focal and reference species
for SP in galGal hetAtr netAur oxyJam stiNae;
do
  cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$SP/$SP.fasta /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$SP/$SP.fasta.fai .
done

# get galGal cDNA fasta 
wget ftp://ftp.ensembl.org/pub/release-100/fasta/gallus_gallus/cdna/Gallus_gallus.GRCg6a.cdna.all.fa.gz
gunzip Gallus_gallus.GRCg6a.cdna.all.fa.gz

# use GMAP to align cDNA to galGal genome
mkdir galGal_db
conda activate gmap 
# build database
gmap_build -D galGal_db -d galGal galGal.fasta
# map cDNA to genome database
gmap -D galGal_db -d galGal -f sampe --sam-extended-cigar -A Gallus_gallus.GRCg6a.cdna.all.fa > galGal.sam

# convert output to a sorted & indexed BAM 
samtools sort galGal.sam -O BAM -o galGal.bam
#| samtools index -b - galGal.sorted.bai
