#!/bin/bash
#SBATCH -J covBeds
#SBATCH -o logs/slurm-%j
#SBATCH -e logs/slurm-%j
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=8000

# run from /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/gatherVCFs_dir/coverage
# sbatch write_coverage_beds.sh spp_name

set -o errexit

cd $1/
mean=$(awk '{sum = sum+$4}{size=size+$2}{avg=sum/size}END{print avg}' $1.summary.tab)

gzip -dc $1.merge.bg.gz | awk -v avg="$mean" -v spp=$1 -f ../sum_cov.awk
