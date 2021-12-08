## in /n/holyscratch01/informatics/swuitchik/ducks/compGen

# copy over orthogroup-specific FASTAs
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/og_fastas .

# align FASTAs
mkdir aligned
sbatch og_fastas/run_mafft.sh

# clean alignments 
sbatch aligned/run_hmmcleaner.sh
