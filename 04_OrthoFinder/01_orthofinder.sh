#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos
# NB: input data is all generated in 03b_CESAR/03_BUSCO.sh

#conda create -c bioconda -n ortho orthofinder diamond

# issues with duplicate sequences in hetAtr; janky grep to remove
cd from_cesar/input_data/
grep -v -f dup1.txt hetAtr.translated.fa | grep -v -f dup2.txt | grep -v -f dup3.txt > hetAtr.translated.dedup.fa

# fix annoying scaffold version issues 
for SP in netAur oxyJam stiNae;
do
	sed '/^>/s/\./_/g' $SP.translated.fa > $SP.translated.clean.fa
done

sed '/^>/s/\./_/g' hetAtr.translated.dedup.fa > hetAtr.translated.clean.fa
cd ../..

# run OrthoFinder 
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/*.clean.fa from_ncbi/input_data/* run_ortho/input_data
sbatch ortho_run.sh
