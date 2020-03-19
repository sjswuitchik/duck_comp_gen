# 4d sites
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/Comparative-Annotation-Toolkit/4d_sites

module load Anaconda3/5.0.1-fasrc02 
source activate py27
module load cufflinks/2.2.1-fasrc01
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread
cd gffread
make release

# get galGal6 annotation from NCBI
cd ..
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# make a CDS-only GFF
gffread/gffread GCF_000002315.6_GRCg6a_genomic.gff --no-pseudo -C -T -o galGal6.filt1.gtf
grep "protein_coding\|CDS" galGal6.filt1.gtf > galGal6.final.gtf

# convert to GP & BED
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod +x ./gtfToGenePred
chmod +x ./genePredToBed

./gtfToGenePred galGal6.final.gtf galGal6.gp
./genePredToBed galGal6.gp galGal6.cds.bed

#run phylo4d from HALtools
module load python/3.6.3-fasrc02
wget http://www.netlib.org/clapack/clapack.tgz
tar -xvzf clapack.tgz
rm clapack.tgz 
mv CLAPACK-3.2.1 clapack
cd clapack
cp make.inc.example make.inc && make f2clib && make blaslib && make lib
export CLAPACKPATH=`pwd`
cd ..
git clone https://github.com/CshlSiepelLab/phast.git
cd phast
cd src && make
cd ../..
export ENABLE_PHYLOP=1

# tree top 1
phast/bin/phyloFit --tree "(((anaPla,(braCan,(ansInd,(ansCyg,ansBra)))),(netAur,(oxyJam,(hetAtr,stiNae)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))))" --subst-mod SSREV --log gallo_top1.out --msa-format MAF ../../CESAR2.0/galloanserae_rooted.maf 

# tree top 2
phast/bin/phyloFit --tree "(((anaPla,(netAur,(oxyJam,(hetAtr,stiNae)))),(braCan,(ansInd,(ansCyg,ansBra)))),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))))" --subst-mod SSREV --log gallo_top2.out --msa-format MAF --out-root galloTop2 ../../CESAR2.0/galloanserae_rooted.maf 

# tree top 3
phast/bin/phyloFit --tree "((((braCan,(ansInd,(ansCyg,ansBra))),(netAur,(oxyJam,(hetAtr,stiNae)))),anaPla),(numMel,(colVir,((tymCupPin,syrMik),(galGal,cotJap)))))" --subst-mod SSREV --log gallo_top3.out --msa-format MAF --out-root galloTop3 ../../CESAR2.0/galloanserae_rooted.maf 

