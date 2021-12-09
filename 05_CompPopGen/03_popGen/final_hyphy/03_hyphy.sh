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

# make BUSTED job scripts
mkdir -p job_scripts_busted/logs

while IFS= read -r file
do
  echo -e '#!/bin/bash' >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -o logs/%j.out" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -e logs/%j.err" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -p shared" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -n 1" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH --mem=9000" >> job_scripts_busted/run_${file}.sh
  echo -e "source activate align\n" >> job_scripts_busted/run_${file}.sh
  echo -e "hyphy busted --alignment ../aligned/${file}.clean.fa --tree ../gene_trees/${file}_tree.txt" >> job_scripts_busted/run_${file}.sh
done < "clean_ogs.tsv"

# create BUSTED batches
cd job_scripts_busted
ls *.sh > scripts
split --numeric-suffixes=1 -l 4000 scripts batch

while IFS= read -r file
do
  sbatch $file
done < "batch01"

while IFS= read -r file
do
  sbatch $file
done < "batch02"

while IFS= read -r file
do
  sbatch $file
done < "batch03"

# make aBSREL job scripts
mkdir -p job_scripts_absrel/logs

while IFS= read -r file
do
  echo -e '#!/bin/bash' >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH -o logs/%j.out" >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH -e logs/%j.err" >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH -p shared" >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH -n 1" >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> job_scripts_absrel/run_${file}.sh
  echo -e "#SBATCH --mem=9000" >> job_scripts_absrel/run_${file}.sh
  echo -e "source activate align\n" >> job_scripts_absrel/run_${file}.sh
  echo -e "hyphy absrel --alignment ../aligned/${file}.clean.fa --tree ../gene_trees/${file}_tree.txt" >> job_scripts_absrel/run_${file}.sh
done < "clean_ogs.tsv"

# create aBSREL batches
cd job_scripts_absrel
ls *.sh > scripts
split --numeric-suffixes=4 -l 4000 scripts batch

while IFS= read -r file
do
  sbatch $file
done < "batch04"

while IFS= read -r file
do
  sbatch $file
done < "batch05"

while IFS= read -r file
do
  sbatch $file
done < "batch06"
