#!/bin/bash
#SBATCH -p holy-info
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 6-00:00:00

set -o nounset -o errexit -o xtrace

export SINGULARITYENV_PYTHONNOUSERSITE=1

readonly SINGULARITY_IMAGE=/n/singularity_images/informatics/cat/cat_v2.1.0-cf7a06e.sif
if [ ! -s cat_work.img ]
then
  readonly tmpdir=$(mktemp -d)
  mkdir -m 777 -p ${tmpdir}/upper
  truncate -s 1T cat_work.img
  singularity exec --cleanenv ${SINGULARITY_IMAGE} mkfs.ext3 -d "${tmpdir}" cat_work.img
  singularity exec --cleanenv --overlay cat_work.img ${SINGULARITY_IMAGE} mkdir /cat_work
  ln -sf /cat_work cat_work
  rm -rf "${tmpdir}"
fi
srun -n 1 env time -v singularity exec \
                           --cleanenv \
                           --no-home \
                           --overlay cat_work.img \
                           "${SINGULARITY_IMAGE}" \
luigi \
  --module cat RunCat \
  --hal=testAll_input/galloanserae.hal \
  --ref-genome=galGal \
  --workers=$((SLURM_JOB_CPUS_PER_NODE-1)) \
  --config=testAll_input/testAll_cat.config \
  --target-genomes='("netAur", "oxyJam", "hetAtr", "stiNae")' \
  --local-scheduler \
  --binary-mode local \
  --augustus \
  --augustus-species=chicken \
  --augustus-cgp \
  --maf-chunksize 850000 \
  --maf-overlap 250000 \
  --assembly-hub
