#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos

conda create -c bioconda -n ortho orthofinder gffread
source activate ortho
mkdir -p from_cesar/input_data
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_* input_data/
cd from_cesar/
# NB: this step creates an an index file for the reference genome in the 2bitdir - if you don't want that to happen, copy the fastas, then run in the input_data dir
gffread cleaned_reordered_hetAtr.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/hetAtr/hetAtr.fasta -y hetAtr.translated.fa -S 
gffread cleaned_reordered_netAur.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/netAur/netAur.fasta -y netAur.translated.fa -S
gffread cleaned_reordered_oxyJam.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/oxyJam/oxyJam.fasta -y oxyJam.translated.fa -S
gffread cleaned_reordered_stiNae.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/stiNae/stiNae.fasta -y stiNae.translated.fa -S
mv *.translated.fa input_data/

# grab example data to check it out
wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
tar zxvf OrthoFinder.tar.gz 
rm OrthoFinder.tar.gz
# make sure install works 
OrthoFinder/orthofinder -f ExampleData/
# tidy
rm -r OrthoFinder/

# run OrthoFinder 
orthofinder -f from_cesar/input_data
