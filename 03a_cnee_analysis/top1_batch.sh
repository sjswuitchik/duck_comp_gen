#!/bin/bash
#SBATCH --mem 60000
#SBATCH -p serial_requeue
#SBATCH -o phylo_top1.out
#SBATCH -e phylo_top1.err
#SBATCH -J phyloAcc
#SBATCH -t 12:00:00

module purge
source setupPhyloAcc.sh
./PhyloAcc gallo_top1.txt
