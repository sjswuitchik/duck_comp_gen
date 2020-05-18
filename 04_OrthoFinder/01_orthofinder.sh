#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos

#conda create -c bioconda -n ortho orthofinder gffread
conda activate ortho

# copy GTFs from CESAR 
mkdir -p from_cesar/input_data
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_* from_cesar/
cd from_cesar/
# translate to protein seq
# NB: this step creates an an index file for the reference genome in the 2bitdir - if you don't want that, copy the fastas over, then run in your working dir
for file in hetAtr netAur oxyJam stiNae;
do
  gffread cleaned_reordered_$file.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta -y $file.translated.fa -S 
  mv $file.translated.fa input_data/
done

# get the rest of the species annotations from the WGA 
#NB not all genomes from WGA are annotated so had to exclude ansBra andInd braCan colVir syrMik and tymCupPin
cd ../
mkdir -p from_ncbi/input_data
cd from_ncbi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/850/225/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/971/095/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/078/875/GCF_002078875.1_NumMel1.0/GCF_002078875.1_NumMel1.0_genomic.gff.gz

gunzip *.gz

# rename the GFFs 
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff ansCyg.gff
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff cotJap.gff
mv GCF_002078875.1_NumMel1.0_genomic.gff numMel.gff
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff anaPla.gff
mv GCF_000002315.6_GRCg6a_genomic.gff galGal.gff

# translate to protein seq
for file in anaPla ansCyg cotJap numMel;
do
  gffread $file.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta -y $file.translated.fa -S
done

# need to strip version information from galGal GFF
python3 stripGFF.py galGal.gff > galGal.stripped.gff
gffread galGal.stripped.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/galGal/galGal.fasta -y galGal.translated.fa -S

# tidy
mv *.translated.fa input_data/

# run OrthoFinder 
cd ..
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/* from_ncbi/input_data/* run_ortho/input_data
orthofinder -f run_ortho/input_data
