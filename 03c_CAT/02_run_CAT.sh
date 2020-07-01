#!/bin/bash
#SBATCH -p holy-info
#SBATCH -N 1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 6-00:00:00
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

singularity exec --cleanenv --home=/var/empty /n/singularity_images/informatics/cat/cat:20200604.sif \
luigi \
  --module cat RunCat \
  --hal=testAll_input/galloanserae.hal \
  --ref-genome=galGal \
  --workers=24 \
  --config=testAll_input/testAll_cat.config \
  --target-genomes='("braCan", "ansInd", "ansCyg", "ansBra", "netAur", "oxyJam", "hetAtr", "stiNae", "numMel", "colVir", "tymCupPin", "syrMik", "cotJap")' \
  --local-scheduler \
  --binary-mode local \
  --augustus \
  --augustus-species=chicken \
  --augustus-cgp \
  --maf-chunksize 850000 \
  --maf-overlap 250000 \
  --assembly-hub
