############################
## Prep & run OrthoFinder ##
############################

# in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021

#conda create -c bioconda conda-forge -n ortho orthofinder diamond perl gffread

mkdir -p run_ortho/input_data
mkdir -p gffs/

# species available on NCBI: anaPla, ansCyg, colVir, cotJap, galGal, numMel downloaded with NCBI Datasets for relevant assemblies

# get duck protein fastas
for file in hetAtr stiNae oxyJam netAur;
do
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/busco/input_data/$file.translated.fa run_ortho/input_data/
done

# get GFFs from Comp Aug annotations for species not available on NCBI & translated to protein FASTA 
for file in braCan tymCupPin syrMik ansInd ansBra;
do
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/$file.gff gffs/
  gffread gffs/$file.gff -g /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/genomes/$file.ncbi.fasta -y run_ortho/input_data/$file.translated.fa -S 
done 

# run OrthoFinder with all species from WGA
sbatch run_orthoFinder.sh

#### results in 
