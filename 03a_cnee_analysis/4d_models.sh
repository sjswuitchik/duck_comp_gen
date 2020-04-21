# 4d sites
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/Comparative-Annotation-Toolkit/4d_sites

module load Anaconda3/2019.10 
source activate py27
module load cufflinks/2.2.1-fasrc01 perl/5.26.1-fasrc01
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
cd gffread
make release

# get galGal6 annotation from NCBI
cd ..
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# clean up GFF to remove partial=true
grep "partial=true" GCF_000002315.6_GRCg6a_genomic.gff | grep "[[:space:]]gene[[:space:]]" | perl -p -e 's/.*(GeneID:\d+;).*/$1/' > geneIDs_remove_parts.txt
grep "gene_biotype=pseudogene" GCF_000002315.6_GRCg6a_genomic.gff | perl -p -e 's/.*(GeneID:\d+;).*/$1/' > geneIDs_remove_pseudo.txt
grep -v -f geneIDs_remove_parts.txt GCF_000002315.6_GRCg6a_genomic.gff | grep -v -f geneIDs_remove_pseudo.txt | grep -v "partial=true" > galGal6.filt.gff 
cp /n/holyscratch01/informatics/swuitchik/ducks_project/cactus_ncbi/galGal.defline.fasta .
gffread/gffread galGal6.filt.gff -J -r -R --no-pseudo -g galGal.defline.fasta -C -o galGal6.filt2.gff
# delete everything in header except ##gff-version 3 b/c for some reason, gffread headers cause errors in gff3ToGenePred

# convert to GP & BED
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod +x ./gff3ToGenePred
chmod +x ./genePredToBed

./gff3ToGenePred galGal6.filt2.gff galGal6.gp
./genePredToBed galGal6.gp galGal6.cds.bed

# extract 4d sites
singularity shell --bind /usr/bin/split --cleanenv /n/singularity_images/informatics/cactus/cactus:2019.03.01--py27hdbcaa40_1.sif 
export PATH=$PWD/phast/bin:$PATH
python /usr/local/lib/python2.7/site-packages/hal/phyloP/halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "(((anaPla,(braCan,(ansInd,(ansCyg,ansBra)))),(netAur,(oxyJam,(hetAtr,stiNae)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae.hal galGal galGal6.cds.bed galloTop1.mod 2> top1.err

<<<<<<< HEAD
python /usr/local/lib/python2.7/site-packages/hal/phyloP/halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "(((anaPla,(netAur,(oxyJam,(hetAtr,stiNae)))),(braCan,(ansInd,(ansCyg,ansBra)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae.hal galGal galGal6.cds.bed galloTop2.mod 2> top2.err
=======
# tree top 1
phast/bin/phyloFit --tree "(((anaPla,(braCan,(ansInd,(ansCyg,ansBra)))),(netAur,(oxyJam,(hetAtr,stiNae)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))))" --subst-mod SSREV --log gallo_top1.out --msa-format MAF --out-root galloTop1 ../../CESAR2.0/galloanserae_rooted.maf 
>>>>>>> 856e0746c379a7e98faecbc1c661b4be5cc7e6b0

python /usr/local/lib/python2.7/site-packages/hal/phyloP/halPhyloPTrain.py \
--numProc 12 \
--noAncestors \
--substMod SSREV \
--tree "((((braCan,(ansInd,(ansCyg,ansBra))),(netAur,(oxyJam,(hetAtr,stiNae)))),anaPla),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))));" \
--targetGenomes anaPla braCan ansInd ansCyg ansBra netAur oxyJam hetAtr stiNae numMel colVir tymCupPin syrMik galGal cotJap \
--precision HIGH \
/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae.hal galGal galGal6.cds.bed galloTop3.mod 2> top3.err

exit

# rename ancestral nodes for PhyloAcc
# conda create -n phast -c bioconda phast
source activate phast 
tree_doctor --name-ancestors galloTop1.mod > galloTop1.named.mod
tree_doctor --name-ancestors galloTop2.mod > galloTop2.named.mod
tree_doctor --name-ancestors galloTop3.mod > galloTop3.named.mod
 

