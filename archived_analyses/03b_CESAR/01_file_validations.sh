######################
## File validations ## 
######################

## archived in /n/holylfs/LABS/informatics/swuitchik/ducks/

## MAF creation and validation

singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif
hal2maf galloanserae.hal galloanserae.maf --refGenome galGal --noAncestors --noDupes

mafTools/bin/mafValidator.py --maf=galloanserae.maf
# no duplication errors from mafValidator.py
mafTools/bin/mafStats --maf galloanserae.maf > stats_test.out 2> stats_test.err

mafTools/bin/mafFilter --maf galloanserae.maf --include hetAtr > hetfilter.out 2> hetfilter.err
mafTools/bin/mafFilter --maf galloanserae.maf --include netAur > netfilter.out 2> netfilter.err
mafTools/bin/mafFilter --maf galloanserae.maf --include stiNae > stifilter.out 2> stifilter.err
mafTools/bin/mafFilter --maf galloanserae.maf --include oxyJam > oxyfilter.out 2> oxyfilter.err


## HAL validation - skip if no issues with previous CESAR runs

# A selection of errors from initial CESAR runs: 

#Genes
#rna-XM_015285278.2
#rna-NM_001030824.1
#rna-XM_001232275.5
#rna-XM_025154377.1
#rna-NM_001127439.1
#rna-XM_015277243.2
#rna-NM_001305257.1
#rna-NM_001306152.2
#rna-XM_025149758.1
#rna-NM_001127439.1
#rna-NM_001278110.2
#rna-XM_025145177.1

grep "ID=rna-XM_015285278.2\|ID=rna-NM_001030824.1\|ID=rna-XM_001232275.5\|ID=rna-XM_025154377.1\|ID=rna-NM_001127439.1\|ID=rna-XM_015277243.2\|ID=rna-NM_001305257.1\|ID=rna-NM_001306152.2\|ID=rna-XM_025149758.1\|ID=rna-NM_001127439.1\|ID=rna-NM_001278110.2\|ID=rna-XM_025145177.1" GCF_000002315.6_GRCg6a_genomic.gff > genecheck.gff
awk -f gff2bed.awk genecheck.gff > genecheck.bed
cat genecheck.bed | python3 genenames.py > genecheck.genes.bed
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed hetAtr hetAtr.genecheck.bed 2> hetgenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed netAur netAur.genecheck.bed 2> netgenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed oxyJam oxyJam.genecheck.bed 2> oxygenecheck.liftover.log
halLiftover --noDupes galloanserae.hal galGal genecheck.genes.bed stiNae stiNae.genecheck.bed 2> stigenecheck.liftover.log
exit

