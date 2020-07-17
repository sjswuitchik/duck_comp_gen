#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos
# NB: input data is all generated in 03b_CESAR/03_BUSCO.sh

#conda create -c bioconda -n ortho orthofinder diamond 
#conda install -c conda-forge perl
conda activate ortho

# issues with duplicate sequences in ducks, need to remove 
cd from_cesar/input_data/
for SP in hetAtr netAur oxyJam stiNae;
do
	perl /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos/rm_dup_acc.pl $SP.translated.fa
done

# fix annoying scaffold version issues 
for SP in hetAtr netAur oxyJam stiNae;
do
	sed '/^>/s/\./_/g' $SP.translated_uniq.fa > $SP.translated.clean.fa
done

cd ../..

# run OrthoFinder 
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/*.clean.fa from_ncbi/input_data/* run_ortho/input_data
sbatch ortho_run.sh
