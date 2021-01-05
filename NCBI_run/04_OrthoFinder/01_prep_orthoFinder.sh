######################
## Prep OrthoFinder ##
######################

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
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/joined_pred/anaPla.gff .
source activate busco
gffread anaPla.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/anaPla.ncbi.fasta -y run_ortho/input_data/anaPla.translated.fa -S

# run OrthoFinder 
sbatch run_orthoFinder.sh
