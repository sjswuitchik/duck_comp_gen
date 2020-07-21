#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH -e ortho-%j.err
#SBATCH -o ortho-%j.out

module purge
module load Anaconda3/2019.10
source activate ortho
orthofinder -o /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/04_OrthoFinder/run_ortho -f run_ortho/input_data
