#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -p shared
#SBATCH --mem=4000
#SBATCH -e ortho-%j.err
#SBATCH -o ortho-%j.out

module purge
module load Anaconda3/2019.10
source activate ortho
orthofinder -o run_ortho/output -f run_ortho/input_data
