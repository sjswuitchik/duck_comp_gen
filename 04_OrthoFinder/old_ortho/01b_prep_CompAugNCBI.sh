## Re-running OrthoFinder with FASTAs from NCBI rather than the Comp Aug translated FASTAs to avoid gene name conversion issues
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/OrthoFinder_jan2021

module load Anaconda/5.0.1-fasrc02
source activate ortho ## conda create command is in 01_prep)

mkdir -p run_ortho/input_data
mkdir -p run_ortho/gffs

# get protein FASTAs for species for which there aren't any publicly available 
for file in hetAtr stiNae oxyJam netAur colVir tymCupPin syrMik braCan ansInd ansBra;
do
  cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/04_OrthoFinder/run_ortho/input_data/$file.translated.fa run_ortho/input_data
done

# get FASTAs for species from NCBI
cd run_ortho/gffs
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/078/875/GCF_002078875.1_NumMel1.0/GCF_002078875.1_NumMel1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/476/345/GCF_015476345.1_ZJU1.0/GCF_015476345.1_ZJU1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/971/095/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz
gunzip *.gz

mv GCF_000002315.6_GRCg6a_genomic.gff galGal.gff
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff ansCyg.gff
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff cotJap.gff
mv GCF_002078875.1_NumMel1.0_genomic.gff numMel.gff
mv GCF_015476345.1_ZJU1.0_genomic.gff anaPla.gff

for file in galGal ansCyg cotJap numMel anaPla;
do
  gffread $file.gff -g /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/genomes/$file.ncbi.fasta -y ../input_data/$file.translated.fa -S
done

# run OrthoFinder with all spp from WGA using NCBI data where available and Comp Aug data where necessary
sbatch run_orthoFinder_v2.sh


