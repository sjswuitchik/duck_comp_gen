################################
## Gene annotation with CESAR ## 
################################

module load Anaconda3/2019.10

# set up and compile in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest

git clone https://github.com/hillerlab/CESAR2.0.git
cd CESAR2.0
make
cd kent/src
make 
cd ../..
export PATH=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest/CESAR2.0/kent/bin:/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest/CESAR2.0/tools:$PATH
mv cesar tools/
cd ..

# organize data needed for variables

# for input genes
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.gp .
python3 stripGPversion.py galGal6.gp > galGal6_stripped.gp

# for 2bit dirs
mkdir 2bitdir
cd 2bitdir
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees/*.fasta .
# strip version number off scaffolds of ref seq
sed '/^>/s/\.*//g' galGal.fasta > galGal_stripped.fasta
rm galGal.fasta
mv galGal_stripped.fasta galGal.fasta
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
./brename -p ".2bit.chrom.sizes" -r ".chrom.sizes" -R

# sort all these files into spp-specific subdirectories in 2bitdir
mv anaPla.* anaPla
mv ansBra.* ansBra
mv ansCyg.* ansCyg
mv ansInd.* ansInd
mv braCan.* braCan
mv colVir.* colVir
mv cotJap.* cotJap
mv galGal.* galGal
mv hetAtr.* hetAtr
mv netAur.* netAur
mv numMel.* numMel
mv oxyJam.* oxyJam
mv stiNae.* stiNae
mv syrMik.* syrMik
mv tymCupPin.* tymCupPin

# for alignment
cd ..
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae_test.maf .
#conda create -c bioconda -n cesar perl perl-scalar-util-numeric
source activate cesar
python3 mafSpeciesScaffoldOnly.py galloanserae_test.maf > galloanserae_stripped.maf
cd CESAR2.0/tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x ./bedToBigBed
chmod +x ./bedSort
cd ../..
export cesarTools=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest/CESAR2.0/tools
export PATH=$PATH:$cesarTools
CESAR2.0/tools/mafIndex galloanserae_stripped.maf galloanserae.bb -chromSizes=2bitdir/galGal/galGal.chrom.sizes

# define variables 
export inputGenes=galGal6.gp
export reference=galGal
export twoBitDir=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest/2bitdir
export alignment=galloanserae.bb
export querySpecies=hetAtr,netAur,oxyJam,stiNae 
export outputDir=CESARoutput 
export resultsDir=geneAnnotation
export maxMemory=50
export profilePath=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesartest/CESAR2.0

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

# tidy
rm -rf $outputDir

