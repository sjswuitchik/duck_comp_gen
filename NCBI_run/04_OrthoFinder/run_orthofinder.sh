#!/bin/bash
#SBATCH -J ortho
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH -p shared
#SBATCH --mem=8000
#SBATCH -e ortho-%j.err
#SBATCH -o ortho-%j.out

#submit from /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021

source activate ortho
orthofinder -o run_ortho/ -f run_ortho/input_data
