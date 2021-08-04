#!/bin/bash
#SBATCH -J sweep
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 01-00:00:00
#SBATCH --mem=12000 

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/sweepfinder
## this is for the Z, W, and MT where the empirical freq spect doesn't load

source activate sweepfinder

for file in CM021763.1.sweep CM021764.1.sweep CM021836.1.sweep;
do
  SweepFinder2 -lg 1000 $file $file.out
done
