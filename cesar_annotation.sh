## CESAR for duck alignment ## 
## https://github.com/hillerlab/CESAR2.0 ## 

## HAL Tools summaries and MAF conversion 
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus on bioinf01

git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit.git
cp ../ducks_cactus/galloanserae.hal Comparative-Annotation-Toolkit/
cd Comparative-Annotation-Toolkit
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
halValidate duck_data/galloanserae.hal
hal2maf galloanserae.hal galloanserae_rooted.maf --refGenome galGal 
halStats data/galloanserae.hal > data/halStats.out
exit

# set up and compile in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus
git clone https://github.com/hillerlab/CESAR2.0.git
cd CESAR2.0
make
cd kent/src
make 
cd ../..
export PATH=`pwd`/kent/bin:`pwd`/tools:$PATH
mv cesar tools/

# organize data needed for variables
# inputGenes (from 4d_sites)
cp ../Comparative-Annotation-Toolkit/4d_sites/galGal6.gp . 

# alignment MAF (from postcactus)
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ducks_cactus/galloanserae_rooted.maf .

# 2bit dir
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo
mkdir 2bitdir
cd 2bitdir
cp -v ../../../ducks_cactus/*.fasta .
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

# alignment
cd ../tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x ./bedToBigBed
chmod +x ./bedSort
cd ..
export cesarTools=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/tools
export PATH=$PATH:$cesarTools
# replace chr names in galGal.chrom.sizes from NC_006089.5 to NC_006089, etc. 
tools/mafIndex galloanserae_rooted.maf galloanserae_rooted.bb -chromSizes=2bitdir/galGal/galGal.chrom.sizes

# define variables 
export inputGenes=galGal6.gp
export reference=galGal
export twoBitDir=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/2bitdir
export alignment=galloanserae_rooted.bb
export querySpecies=hetAtr,netAur,oxyJam,stiNae 
export outputDir=CESARoutput 
export resultsDir=geneAnnotation
export maxMemory=50
export profilePath=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/

# create CESAR 

#### 
# these are weird rearrangements I had to do to get the jobList to run properly - try without these (ie/ go straight to formatGenePred.pl) and see if it works first. If not, may need these rearrangements
cp tools/cesar extra/
cd extra/tables/human/
chmod +x *.txt
cd ../../
mkdir extra/
mv tables/ extra/ 
####
formatGenePred.pl ${inputGenes} ${inputGenes}.forCESAR ${inputGenes}.discardedTranscripts -longest
for transcript in `cut -f1 ${inputGenes}.forCESAR`; do 
   echo "annotateGenesViaCESAR.pl ${transcript} ${alignment} ${inputGenes}.forCESAR ${reference} ${querySpecies} ${outputDir} ${twoBitDir} ${profilePath} -maxMemory ${maxMemory}"
done > jobList

# realign all genes (this step with test data takes 10.5 hr on single core)
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

