#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 55000
#SBATCH -p shared
#SBATCH -o phylo_top2.out
#SBATCH -e phylo_top2.err
#SBATCH -J phyloAcc
#SBATCH -t 12:00:00

module purge
source setupPhyloAcc.sh
./PhyloAcc gallo_top2.txt
