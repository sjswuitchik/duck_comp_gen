#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p shared
#SBATCH -e compAug_denovo_%A.e
#SBATCH -o compAug_denovo_%A.o
#SBATCH -J compAug_denovo
#SBATCH --mem=64000
#SBATCH -t 23:00:00

module purge
module load Anaconda/5.0.1-fasrc01
source activate compAug

for ali in *.maf;
do
id=${ali%.maf} 
augustus \
--species=chicken \
--softmasking=1 \
--treefile=top1.nwk \
--alnfile=$ali \
--dbaccess=../chicken.db \
--speciesfilenames=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes.tbl \
--alternatives-from-evidence=0 \
--/CompPred/outdir=pred$id > aug$id.out 2> err$id.out &
done
