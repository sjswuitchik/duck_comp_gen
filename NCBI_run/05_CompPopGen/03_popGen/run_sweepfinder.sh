#!/bin/bash
#SBATCH -J sweep
#SBATCH -o out
#SBATCH -e err
#SBATCH -p bigmem
#SBATCH -n 1
#SBATCH -t 12-00:00:00
#SBATCH --mem=190G

# submit from /n/holyscratch01/informatics/swuitchik/ducks/compGen/sweepfinder

source activate sweepfinder

for file in CM*.sweep;
do
  SweepFinder2 -lg 1000 $file hetAtr.spect $file.out
done
