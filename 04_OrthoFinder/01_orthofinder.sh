#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos
# NB: input data is all generated in 03b_CESAR/03_BUSCO.sh

#conda create -c bioconda -n ortho orthofinder diamond

# issues with duplicate sequences in hetAtr; use seqtk to remove seq
cd from_cesar/input_data
git clone https://github.com/lh3/seqtk.git;
cd seqtk
make
cd ..
seqtk/seqtk subseq hetAtr.translated.fa dups.list > hetAtr.translated.dedups.fa

# fix annoying scaffold version issues 
cd from_cesar/input_data/
for SP in netAur oxyJam stiNae;
do
	sed '/^>/s/\./_/g' $SP.translated.fa > $SP.translated.clean.fa
done





cd ../..

# run OrthoFinder 
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/*.clean.fa from_ncbi/input_data/* run_ortho/input_data
sbatch ortho_run.sh
