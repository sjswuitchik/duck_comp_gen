######################
## Prep OrthoFinder ##
######################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/04_OrthoFinder

#conda create -c bioconda -n ortho orthofinder diamond 
#conda install -c conda-forge perl
conda activate ortho

mkdir -p run_ortho/input_data
cd run_ortho/input_data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/077/185/GCF_011077185.1_BPBGC_Ojam_1.0/GCF_011077185.1_BPBGC_Ojam_1.0_protein.faa.gz
gunzip *.gz

# run OrthoFinder 
sbatch ortho_run.sh
