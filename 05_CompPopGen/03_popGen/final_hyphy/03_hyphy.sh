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

# remove duplicate sequences and produce combined alignment + tree for HYPHY input
cd ..
mkdir clean_aligned
sbatch run_removeDups.sh

# create list of OGs that were dedup'd
cd clean_aligned 
ls *.nh > uniq
sed -i 's/\_uniq\.nh//g' uniq
mv uniq ..
Rscript split_dups.R

## make BUSTED job scripts
mkdir -p job_scripts_busted/logs

# make scripts for runs without duplicate seqs
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
  echo -e "hyphy busted --alignment aligned/${file}_nuc.fa_hmm.fasta.filtered --tree gene_trees/${file}_tree.txt" >> job_scripts_busted/run_${file}.sh
done < "split"

# make scripts for runs with duplicate seqs 
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
  echo -e "hyphy busted --alignment clean_aligned/${file}_uniq.nh" >> job_scripts_busted/run_${file}.sh
done < "uniq"

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

## make aBSREL job scripts
cd ..
mkdir -p job_scripts_absrel/logs

# make scripts for runs without duplicate seqs
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
  echo -e "hyphy absrel --alignment aligned/${file}_nuc.fa_hmm.fasta.filtered --tree gene_trees/${file}_tree.txt" >> job_scripts_absrel/run_${file}.sh
done < "split"

# make scripts for runs with duplicate seqs 
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
  echo -e "hyphy absrel --alignment clean_aligned/${file}_uniq.nh" >> job_scripts_absrel/run_${file}.sh
done < "uniq"

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
