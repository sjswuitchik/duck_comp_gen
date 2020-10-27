#!/usr/bin/bash

#SBATCH -t 800
#SBATCH --mem 40000
#SBATCH -p serial_requeue,bos-info
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --array=1-50

DATAPATH=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc/go_perms
Rscript cnee_go_perms.R ${SLURM_ARRAY_TASK_ID} 5 50 $DATAPATH 
