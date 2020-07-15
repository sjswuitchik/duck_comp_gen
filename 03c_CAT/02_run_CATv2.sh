#!/bin/sh
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --partition=shared
#SBATCH --time=06-00:00:00
set -o nounset -o errexit -o xtrace
readonly SINGULARITY_IMAGE=/n/singularity_images/informatics/cat/cat_v2.1.0-cf7a06e.sif

mkdir -p cat_work /scratch/cat-toil.${SLURM_JOB_ID}
ln -sf /scratch/cat-toil.${SLURM_JOB_ID} cat_work/toil

# this should probably be set in the CAT image
export SINGULARITYENV_PYTHONNOUSERSITE=1

srun -n 1 env time -v singularity exec \
                           --cleanenv \
                           --no-home \
                           "${SINGULARITY_IMAGE}" \
  luigi \
    --module cat RunCat \
    --hal=testAll_input/galloanserae.hal \
    --ref-genome=galGal \
    --workers=$((SLURM_JOB_CPUS_PER_NODE-1)) \
    --config=testAll_input/testAll_cat.config \
    --local-scheduler \
    --binary-mode local \
    --augustus \
    --augustus-species=chicken \
    --augustus-cgp \
    --assembly-hub 

rm -rf /scratch/cat-toil.${SLURM_JOB_ID}
