################################
## Gene annotation with CESAR ## 
################################

module load Anaconda3/2019.10
#conda create -c bioconda -n cesar perl perl-scalar-util-numeric
source activate cesar

# set up and compile in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar

git clone https://github.com/hillerlab/CESAR2.0.git
cd CESAR2.0
make
cd kent/src
make 
cd ../..
export PATH=`pwd`/kent/bin:`pwd`/tools:$PATH
export profilePath=`pwd`
cp cesar tools/
cd ..

# organize data needed for variables

# for input genes
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.gp .
# strip version number off scaffolds in GenePred
python3 stripGPversion.py galGal6.gp > galGal6_stripped.gp

# for 2bit dirs
mkdir 2bitdir
cd 2bitdir
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees/*.fasta .
# strip version number off scaffolds in ref seq
awk -f stripFasta.awk galGal.defline.fasta > galGal_stripped.fasta
mv galGal_stripped.fasta galGal.fasta
rm galGal.defline.fasta

# these loops are ugly but quick - adjust file names with brename after
for file in *.fasta; 
do
	./faToTwoBit $file $file.2bit
done

for file in *.2bit;
do
	./twoBitInfo $file stdout | sort -k2rn > $file.chrom.sizes
done

for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	mkdir $file
done

wget https://github.com/shenwei356/brename/releases/download/v2.10.0/brename_linux_amd64.tar.gz
tar zxvf brename_linux_amd64.tar.gz 
rm brename_linux_amd64.tar.gz 
chomd +x ./brename
./brename -p ".defline." -r "." -R
./brename -p ".fasta.2bit" -r ".2bit" -R

# sort all these files into spp-specific subdirectories in 2bitdir
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	mv $file.* $file/
done

# each '*.2bit.chrom.sizes' file needs to be named 'chrom.sizes', so need to rename after sorting 
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	mv $file/$file.2bit.chrom.sizes $file/chrom.sizes
done

# for alignment
cd ..
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae_test.maf .
python3 mafSpeciesScaffoldOnly.py galloanserae_test.maf > galloanserae_stripped.maf
cd CESAR2.0/tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x ./bedToBigBed
chmod +x ./bedSort
cd ..
export cesarTools=`pwd`/tools
export PATH=$PATH:$cesarTools
cd ..
CESAR2.0/tools/mafIndex galloanserae_stripped.maf galloanserae.bb -chromSizes=2bitdir/galGal/chrom.sizes

# define variables 
export inputGenes=galGal6_stripped.gp
export reference=galGal
export twoBitDir=2bitdir
export alignment=galloanserae.bb
export querySpecies=hetAtr,netAur,oxyJam,stiNae 
export outputDir=CESARoutput 
export resultsDir=geneAnnotation
export maxMemory=50

chmod +x CESAR2.0/tools/formatGenePred.pl
formatGenePred.pl ${inputGenes} ${inputGenes}.forCESAR ${inputGenes}.discardedTranscripts -longest
for transcript in `cut -f1 ${inputGenes}.forCESAR`; do 
   echo "annotateGenesViaCESAR.pl ${transcript} ${alignment} ${inputGenes}.forCESAR ${reference} ${querySpecies} ${outputDir} ${twoBitDir} ${profilePath} -maxMemory ${maxMemory}"
done > jobList

# realign all genes
chmod +x jobList
./jobList > jobList.out 2> jobList.err

# convert results
for species in `echo $querySpecies | sed 's/,/ /g'`; do 
  echo "bed2GenePred.pl $species $outputDir /dev/stdout | awk '{if (\$4 != \$5) print \$0}' > $resultsDir/$species.gp"
done > jobListGenePred
chmod +x jobListGenePred
mkdir $resultsDir
./jobListGenePred

# tidy up
rm -rf $outputDir

# convert each query species output to GTF for filtering
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod +x ./genePredToGtf

#genePredToGtf file refGene.input hg19refGene.gtf
genePredToGtf geneAnnotation/hetAtr.gp ____ hetAtr.gtf
