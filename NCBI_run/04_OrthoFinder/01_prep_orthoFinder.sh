######################
## Prep OrthoFinder ##
######################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/04a_OrthoFinder_compAug

#conda create -c bioconda -n ortho orthofinder diamond 
#conda install -c conda-forge perl
conda activate ortho

mkdir -p run_ortho/input_data
for file in hetAtr stiNae oxyJam netAur;
do
  cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/busco/input_data/$file.translated.fa run_ortho/input_data/
done

# run OrthoFinder 
sbatch 02_ncbi_ortho.sh
