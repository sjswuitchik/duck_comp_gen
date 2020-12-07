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

cd maflinks/

for ali in $(seq 1 962);
do
 id=$ali # remove .maf suffix
 augustus \
  --species=chicken \
  --softmasking=1 \
  --treefile=top1.nwk \
  --alnfile=$ali \
  --dbaccess=chicken_rnaseq.db \
  --speciesfilenames=../genomes.tbl \
  --alternatives-from-evidence=0 \
  --dbhints=1 \
  --UTR=1 \
  --allow_hinted_splicesites=atac \
  --extrinsicCfgFile=extrinsic-rnaseq.cfg \
  --/CompPred/outdir=pred$id > aug$id.out 2> err$id.out &
done
