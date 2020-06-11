##############################
## Gene annotation with CAT ##
##############################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/
# also testing in /scratch/swuitchik on bioinf01 for optimization and debugging
# https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit

module load samtools/1.10-fasrc01

# set up CAT
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit.git
pip install -e Comparative-Annotation-Toolkit

# extract only focal and reference species from WGA
mkdir -p Comparative-Annotation-Toolkit/input_data
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
hal2maf galloanserae.hal galloForCAT.maf --refGenome galGal --noAncestors --noDupes --targetGenomes galGal,hetAtr,netAur,oxyJam,stiNae
maf2hal galloForCAT.maf galloForCAT.hal --refGenome galGal
exit 

# set up input data
mv galloForCAT.hal Comparative-Annotation-Toolkit/input_data
cd Comparative-Annotation-Toolkit/input_data/
# get filtered GFF for galGal
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.filtpy.gff . 

# get galGal protein alignment
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_protein.faa.gz

# get indexed fastas for focal and reference species
for SP in galGal hetAtr netAur oxyJam stiNae;
do
  cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$SP/$SP.fasta /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$SP/$SP.fasta.fai .
done

cd ..

# run CAT
singularity exec --cleanenv --no-home /n/singularity_images/informatics/cat/cat:20200604.sif \
luigi \
  --module cat RunCat \
  --hal=input_data/galloForCAT.hal \
  --ref-genome=galGal \
  --workers=10 \
  --config=input_data/cat.config \
  --local-scheduler \
  --binary-mode local \
  --augustus \
  --augustus-cgp \
  --augustus-pb \
  --assembly-hub
