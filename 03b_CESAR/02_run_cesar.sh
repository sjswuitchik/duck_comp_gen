################################
## Gene annotation with CESAR ## 
################################

module load Anaconda3/2019.10

# set up and compile in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus

git clone https://github.com/hillerlab/CESAR2.0.git
cd CESAR2.0
make
cd kent/src
make 
cd ../..
export PATH=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/kent/bin:/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/tools:$PATH
mv cesar tools/
cd ..

# organize data needed for variables

# for input genes
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.gp .

# for 2bit dirs
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo
mkdir 2bitdir
cd 2bitdir
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees/*.fasta .
# these loops are ugly but quick - adjust file names with brename after
for file in *.fasta; 
do
	faToTwoBit $file $file.2bit
done

for file in *.2bit;
do
	twoBitInfo $file stdout | sort -k2rn > $file.chrom.sizes
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
./brename -p "oxyJam.masked" -r "oxyJam" -R
./brename -p "stiNae.masked" -r "stiNae" -R
./brename -p "hetAtr.masked" -r "hetAtr" -R

# sort all these files into spp-specific subdirectories in 2bitdir

# for alignment
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/galloanserae.maf .
#conda create -c bioconda -n cesar perl perl-scalar-util-numeric
source activate cesar
python3 mafSpeciesScaffoldOnly.py galloanserae.maf > galloanserae_stripped.maf
cd ../tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x ./bedToBigBed
chmod +x ./bedSort
cd ..
export cesarTools=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/tools
export PATH=$PATH:$cesarTools
tools/mafIndex galloanserae_stripped.maf galloanserae.bb -chromSizes=2bitdir/galGal/galGal.chrom.sizes

#### 
# these are weird rearrangements I had to do to get the jobList to run properly - try without these (ie/ go straight to formatGenePred.pl) and see if it works first. If not, may need these rearrangements
cp tools/cesar extra/
cp tools/cesar .
cd extra/tables/human/
chmod +x *.txt
cd ../../
####

# define variables 
export inputGenes=galGal6.gp
export reference=galGal
export twoBitDir=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/2bitdir
export alignment=galloanserae.bb
export querySpecies=hetAtr,netAur,oxyJam,stiNae 
export outputDir=CESARoutput 
export resultsDir=geneAnnotation
export maxMemory=50
export profilePath=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0

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

