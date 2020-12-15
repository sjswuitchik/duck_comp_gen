#!/bin/sh
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --partition=holy-info
#SBATCH --time=14-00:00:00
set -o nounset -o errexit -o xtrace

readonly SINGULARITY_IMAGE=/n/singularity_images/informatics/augustus/augustus_3.4.0-afreedman-build.sif

srun -n 1 env time -v singularity exec \
                           --cleanenv \
                           --no-home \
                           "${SINGULARITY_IMAGE}" \
for ali in *.maf;
do
 augustus \
  --species=chicken \
  --softmasking=1 \
  --treefile=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/top1.nwk \
  --alnfile=$ali \
  --dbaccess=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/chicken_rnaseq.db \
  --speciesfilenames=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes.tbl \
  --alternatives-from-evidence=0 \
  --dbhints=1 \
  --extrinsicCfgFile=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/extrinsic-rnaseq.cfg \
  --/CompPred/outdir=pred$id
done


