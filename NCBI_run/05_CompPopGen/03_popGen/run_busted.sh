#!/bin/bash
#SBATCH -J busted
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 06-23:00:00
#SBATCH --mem=12000

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/

source activate busted

# hyphy busted --help for man page
