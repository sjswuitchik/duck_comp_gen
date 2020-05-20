#!/bin/sh
#SBATCH --nodes=1
# allow use of all the memory on the node
#SBATCH --mem=0
# request all CPU cores on the node
#SBATCH --exclusive
# Customize --time --partition as appropriate
#SBATCH --time=14-00:00:00
#SBATCH --partition=holy-info


set -o nounset -o errexit -o xtrace

########################################
# parameters
########################################
readonly CACTUS_IMAGE=/n/singularity_images/informatics/cactus/cactus:2019-11-29.sif
readonly JOBSTORE_IMAGE=jobStore.img # cactus jobStore; will be created if it doesn't exist
readonly SEQFILE=gallo_ncbi.txt
readonly OUTPUTHAL=gallo_ncbi.hal

########################################
# ... don't modify below here ...

readonly CACTUS_SCRATCH=/scratch/cactus-${SLURM_JOB_ID}

if [ ! -e "${JOBSTORE_IMAGE}" ]
then
  restart=''
  mkdir -m 777 -p ${CACTUS_SCRATCH}/upper
  truncate -s 1T "${JOBSTORE_IMAGE}"
  singularity exec /n/singularity_images/informatics/cactus/e2fsprogs:1.45.2.sif mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"
else
  restart='--restart'
fi

# Use empty home & /tmp directories in the container (to avoid, e.g., pip-installed packages in ~/.local)
mkdir -m 700 -p ${CACTUS_SCRATCH}/home ${CACTUS_SCRATCH}/tmp

# the toil workDir must be on the same file system as the cactus jobStore
singularity exec --overlay ${JOBSTORE_IMAGE} ${CACTUS_IMAGE} mkdir -p /cactus/workDir
srun -n 1 singularity exec --cleanenv \
                           --overlay ${JOBSTORE_IMAGE} \
                           --bind ${CACTUS_SCRATCH}/tmp:/tmp \
                           --home ${CACTUS_SCRATCH}/home:/home/cactus ${CACTUS_IMAGE} \
  cactus ${CACTUS_OPTIONS-} ${restart-} --workDir=/cactus/workDir --binariesMode local /cactus/jobStore "${SEQFILE}" "${OUTPUTHAL}"

# /tmp would eventually be purged, but just in case the
# next job to run on this node needs lots of /space...

rm -rf ${CACTUS_SCRATCH}
