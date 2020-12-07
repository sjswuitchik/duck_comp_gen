#!/bin/bash
#SBATCH -p shared
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --mem 32000
#SBATCH --array 0-187

source /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc_control/PhyloAcc/setupPhyloAcc.sh
./PhyloAcc top3_param/run${SLURM_ARRAY_TASK_ID}