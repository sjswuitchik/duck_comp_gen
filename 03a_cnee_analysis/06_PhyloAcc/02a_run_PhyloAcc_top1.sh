#!/bin/bash
#SBATCH -p shared
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --mem 32000
#SBATCH --array 0-178

source /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/PhyloAcc/setupPhyloAcc.sh
./PhyloAcc top1_param/run${SLURM_ARRAY_TASK_ID}
