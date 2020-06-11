##############################
## Create 4d neutral models ## 
##############################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis

module load Anaconda3/2019.10 

mkdir 4d_neut_models
cd 4d_neut_models/

# get galGal6 annotation from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# clean up GFF to remove partial=true
grep "partial=true" GCF_000002315.6_GRCg6a_genomic.gff | grep "[[:space:]]gene[[:space:]]" | perl -p -e 's/.*(GeneID:\d+).*/$1/' > geneIDs_remove_parts.txt
grep "gene_biotype=protein_coding" GCF_000002315.6_GRCg6a_genomic.gff | grep "[[:space:]]gene[[:space:]]" | perl -p -e 's/.*(GeneID:\d+).*/$1/' > geneIDs_keep_prot.txt
grep -v -f geneIDs_remove_parts.txt GCF_000002315.6_GRCg6a_genomic.gff | grep -f geneIDs_keep_prot.txt > galGal6.filt.gff
sed -i '1i#!gff-spec-version 1.21' galGal6.filt.gff && sed -i '1i##gff-version 3' galGal6.filt.gff
python3 CustomExtractPassFiltMrnaFromGff.py
python3 WriteFilteredGff.py

# convert to GP & BED (for use in CESAR as well) 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod +x ./gff3ToGenePred
chmod +x ./genePredToBed

./gff3ToGenePred galGal6.filtpy.gff galGal6.gp
./genePredToBed galGal6.gp galGal6.cds.bed

# extract 4d sites
singularity exec --cleanenv /n/singularity_images/informatics/hal/hal:20200527.sif halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "(((anaPla,(braCan,(ansInd,(ansCyg,ansBra)))),(netAur,(oxyJam,(hetAtr,stiNae)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/gallo_ncbi.hal galGal galGal6.cds.bed galloTop1.mod 2> top1.err

singularity exec --cleanenv /n/singularity_images/informatics/hal/hal:20200527.sif halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "(((anaPla,(netAur,(oxyJam,(hetAtr,stiNae)))),(braCan,(ansInd,(ansCyg,ansBra)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/gallo_ncbi.hal galGal galGal6.cds.bed galloTop2.mod 2> top2.err

singularity exec --cleanenv /n/singularity_images/informatics/hal/hal:20200527.sif halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "((((braCan,(ansInd,(ansCyg,ansBra))),(netAur,(oxyJam,(hetAtr,stiNae)))),anaPla),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/gallo_ncbi.hal galGal galGal6.cds.bed galloTop3.mod 2> top3.err

exit

# rename ancestral nodes for PhyloAcc
#conda create -n phast -c bioconda phast
source activate phast 
tree_doctor --name-ancestors galloTop1.mod > galloTop1.named.mod
tree_doctor --name-ancestors galloTop2.mod > galloTop2.named.mod
tree_doctor --name-ancestors galloTop3.mod > galloTop3.named.mod
