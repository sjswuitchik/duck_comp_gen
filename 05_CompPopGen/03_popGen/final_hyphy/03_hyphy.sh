# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

## make BUSTED job scripts
mkdir -p job_scripts_busted/logs


# make scripts for runs without duplicate seqs
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -o logs/${file}.out" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -e logs/${file}.err" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -p shared" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -n 1" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> job_scripts_busted/run_${file}.sh
  echo -e "conda activate align\n" >> job_scripts_busted/run_${file}.sh
  echo -e "hyphy busted --alignment og_fastas/${file}_nuc.fa_codon.msa --tree og_fastas/${file}_tree.txt" >> job_scripts_busted/run_${file}.sh
done < "og_fastas/split_aligns"

# make scripts for runs with duplicate seqs 
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -o logs/${file}.out" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -e logs/${file}.err" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -p shared" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -n 1" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> job_scripts_busted/run_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> job_scripts_busted/run_${file}.sh
  echo -e "conda activate align\n" >> job_scripts_busted/run_${file}.sh
  echo -e "hyphy busted --alignment og_fastas/${file}_uniq.nh" >> job_scripts_busted/run_${file}.sh
done < "og_fastas/uniq"

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
  echo -e "hyphy absrel --alignment og_fastas/${file}_nuc.fa_codon.msa --tree og_fastas/${file}_tree.txt" >> job_scripts_absrel/run_${file}.sh
done < "og_fastas/split"

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
  echo -e "hyphy absrel --alignment og_fastas/${file}_uniq.nh" >> job_scripts_absrel/run_${file}.sh
done < "og_fastas/uniq"

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
