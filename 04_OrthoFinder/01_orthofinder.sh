#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos
# NB: input data is all generated in 03b_CESAR/03_BUSCO.sh

#conda create -c bioconda -n ortho orthofinder
conda activate ortho

# run OrthoFinder 
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/* from_ncbi/input_data/* run_ortho/input_data
orthofinder -o run_ortho/output -f run_ortho/input_data
