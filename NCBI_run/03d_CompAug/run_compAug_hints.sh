#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p shared
#SBATCH -e compAug_hints_%A.e
#SBATCH -o compAug_hints_%A.o
#SBATCH -J compAug_hints
#SBATCH --mem=64000
#SBATCH -t 5-00:00:00

module purge
module load Anaconda/5.0.1-fasrc01
source activate compAug

for ali in *.maf;
do
 id=$ali # remove .maf suffix
 augustus \
  --species=chicken \
  --softmasking=1 \
  --treefile=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/top1.nwk \
  --alnfile=$ali \
  --dbaccess=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/chicken_rnaseq.db \
  --speciesfilenames=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes.tbl \
  --alternatives-from-evidence=0 \
  --dbhints=1 \
  --UTR=1 \
  --allow_hinted_splicesites=atac \
  --extrinsicCfgFile=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/extrinsic-rnaseq.cfg \
  --/CompPred/outdir=pred$id > aug$id.out 2> err$id.out &
done
