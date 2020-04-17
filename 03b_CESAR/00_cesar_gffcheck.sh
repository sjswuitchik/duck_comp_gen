###############################################################################################################
## If you are getting errors from CESAR, here are some steps for troubleshooting the HAL, MAF, and GFF files ## 
###############################################################################################################

# Selection of errors from initial CESAR runs

cesar_52771093_320.out
cesar_52772274_601.out
cesar_52772274_688.out
cesar_52774212_1190.out
cesar_52782332_1516.out
cesar_52902022_2096.out
cesar_52902022_2466.out
cesar_52913403_3404.out
cesar_52918789_3992.out
cesar_52918789_4151.out
cesar_52918789_4230.out
cesar_52918789_4406.out


Genes
rna-XM_015285278.2
rna-NM_001030824.1
rna-XM_001232275.5
rna-XM_025154377.1
rna-NM_001127439.1
rna-XM_015277243.2
rna-NM_001305257.1
rna-NM_001306152.2
rna-XM_025149758.1
rna-NM_001127439.1
rna-NM_001278110.2
rna-XM_025145177.1

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/cesar.gffcheck 

grep "ID=rna-XM_015285278.2\|ID=rna-NM_001030824.1\|ID=rna-XM_001232275.5\|ID=rna-XM_025154377.1\|ID=rna-NM_001127439.1\|ID=rna-XM_015277243.2\|ID=rna-NM_001305257.1\|ID=rna-NM_001306152.2\|ID=rna-XM_025149758.1\|ID=rna-NM_001127439.1\|ID=rna-NM_001278110.2\|ID=rna-XM_025145177.1" GCF_000002315.6_GRCg6a_genomic.gff > genecheck.gff
awk -f gff2bed.awk genecheck.gff > genecheck.bed
cat genecheck.bed | python3 genenames.py > genecheck.genes.bed
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed hetAtr hetAtr.genecheck.bed 2> hetgenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed netAur netAur.genecheck.bed 2> netgenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed oxyJam oxyJam.genecheck.bed 2> oxygenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed stiNae stiNae.genecheck.bed 2> stigenecheck.liftover.log
exit

cp ../galloanserae_final.maf .

# mafTools
source activate py27
# conda install -c conda-forge scipy numpy 
module load gcc/8.2.0-fasrc01 
git clone https://github.com/dentearl/mafTools.git
git clone https://github.com/benedictpaten/sonLib.git
cd sonLib
make
cd ..
git clone https://github.com/benedictpaten/pinchesAndCacti.git
cd pinchesAndCacti
make
cd ../mafTools
make
cd ..
# no errors when run with the --ignore flag
mafTools/bin/mafValidator.py --maf=galloanserae_final.maf --ignoreDuplicateColumns > val.log 2> val.err
# duplicate column error without the --ignore flag
mafTools/bin/mafValidator.py --maf=galloanserae_final.maf > val.log 2> val.err
mafTools/bin/mafStats --maf galloanserae_final.maf > stats.log 2> stats.err
# mafTools make didn't actually work for everything - needed to add `-lm` to the end of line 25 and 30 in the mafExtrator makefile
# made need to `make clean` then `make` again
mafTools/bin/mafDuplicateFilter --maf galloanserae_final.maf > galloanserae_pruned.maf
mafTools/bin/mafValidator.py --maf=galloanserae_pruned.maf > val_pruned.log 2> val_pruned.err 
mafTools/bin/mafDuplicateFilter --maf galloanserae_pruned.maf > galloanserae_pruned2.maf
mafTools/bin/mafValidator.py --maf=galloanserae_pruned2.maf > val_pruned2.log 2> val_pruned2.err 
