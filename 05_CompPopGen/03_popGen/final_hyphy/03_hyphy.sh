# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

# copy over gene trees
mkdir gene_trees
cut -f1 /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/clean_orthos.tsv | sed -e "1d" > clean_ogs.tsv

while IFS= read -r file
do
  cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Gene_Trees/${file}_tree.txt gene_trees/
done < "clean_ogs.tsv"

# edit gene trees to match alignments
cd gene_trees
for file in *.txt;
do
  sed -i 's/\_protein//g' $file 
  sed -i 's/\_translated//g' $file
done

# make job scripts
mkdir job_scripts
mkdir logs

while IFS= read -r file
do
  echo -e '#!/bin/bash' >> job_scripts/run_${file}.sh
  echo -e "#SBATCH -o logs/%j.out" >> job_scripts/run_${file}.sh
  echo -e "#SBATCH -e logs/%j.err" >> job_scripts/run_${file}.sh
  echo -e "#SBATCH -p shared" >> job_scripts/run_${file}.sh
  echo -e "#SBATCH -n 1" >> job_scripts/run_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> job_scripts/run_${file}.sh
  echo -e "#SBATCH --mem=9000" >> job_scripts/run_${file}.sh
  echo -e "source activate align\n" >> job_scripts/run_${file}.sh
  echo -e "hyphy busted --alignment aligned/${file}_nuc.fa_hmm.fasta --tree gene_trees/${file}_tree.txt" >> job_scripts/run_${file}.sh
done < "clean_ogs.tsv"

# create batches
ls job_scripts/*.sh > job_scripts/scripts
split -l 4000 job_scripts/scripts job_scripts/batch

for file in job_scripts/batch*;
do
  sed -i 's/run/sbatch\ run/g' $file
done



