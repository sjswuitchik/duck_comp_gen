#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 55000
#SBATCH -p shared
#SBATCH -o phylo_top3.out
#SBATCH -e phylo_top3.err
#SBATCH -J phyloAcc
#SBATCH -t 12:00:00

module purge
source setupPhyloAcc.sh
./PhyloAcc gallo_top3.txt
