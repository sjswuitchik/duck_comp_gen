############################
## Prep & run OrthoFinder ##
############################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/04a_OrthoFinder_compAug

#conda create -c bioconda -n ortho orthofinder diamond 
#conda install -c conda-forge perl

mkdir -p run_ortho/input_data
# get ducks & requisite outgroup protein fastas
for file in hetAtr stiNae oxyJam netAur galGal;
do
  cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/busco/input_data/$file.translated.fa run_ortho/input_data/
done
# grab mallard GFF and translate to protein fasta
mkdir gffs/
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/joined_pred/anaPla.gff gffs/
source activate busco
gffread gffs/anaPla.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/anaPla.ncbi.fasta -y run_ortho/input_data/anaPla.translated.fa -S

# run OrthoFinder with only two outgroups to test
sbatch run_orthoFinder.sh

### results in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/04a_OrthoFinder_compAug/run_ortho/Results_Jan05

# grab other species from WGA and translate to protein FASTAs
for file in ansBra ansCyg ansInd braCan colVir cotJap numMel syrMik tymCupPin;
do
  cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/joined_pred/$file.gff gffs/
  gffread gffs/$file.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/$file.ncbi.fasta -y run_ortho/input_data/$file.translated.fa -S 
done 

# run OrthoFinder with all species from WGA
sbatch run_orthoFinder.sh

#### results in 
