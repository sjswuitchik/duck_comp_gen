#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 4000
#SBATCH -t 02:00:00
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

singularity exec --cleanenv --home=/var/empty /n/singularity_images/informatics/cat/cat:20200604.sif \
luigi \
  --module cat RunCat \
  --hal=input_data/galloForCAT.hal \
  --ref-genome=galGal \
  --workers=24 \
  --config=input_data/cat.config \
  --target-genomes='("hetAtr", "netAur", "oxyJam", "stiNae")' \
  --local-scheduler \
  --binary-mode local \
  --augustus \
  --augustus-species=chicken \
  --augustus-cgp \
  --maf-chunksize 850000 \
  --maf-overlap 250000 \
  --assembly-hub > stdout.log 2> stdout.err
