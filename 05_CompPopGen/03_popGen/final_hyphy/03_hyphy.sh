# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/og_fastas

## make BUSTED job scripts
mkdir logs

# make scripts for runs without duplicate seqs
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> run_busted_${file}.sh
  echo -e "#SBATCH -o logs/${file}_busted.out" >> run_busted_${file}.sh
  echo -e "#SBATCH -e logs/${file}_busted.err" >> run_busted_${file}.sh
  echo -e "#SBATCH -p shared" >> run_busted_${file}.sh
  echo -e "#SBATCH -n 1" >> run_busted_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> run_busted_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> run_busted_${file}.sh
  echo -e "source activate align\n" >> run_busted_${file}.sh
  echo -e "hyphy busted --alignment ${file}_nuc.fa_codon.msa --tree ${file}_tree.txt" >> run_busted_${file}.sh
done < "split_aligns"

# make scripts for runs with duplicate seqs 
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> run_busted_${file}.sh
  echo -e "#SBATCH -o logs/${file}_busted.out" >> run_busted_${file}.sh
  echo -e "#SBATCH -e logs/${file}_busted.err" >> run_busted_${file}.sh
  echo -e "#SBATCH -p shared" >> run_busted_${file}.sh
  echo -e "#SBATCH -n 1" >> run_busted_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> run_busted_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> run_busted_${file}.sh
  echo -e "source activate align\n" >> run_busted_${file}.sh
  echo -e "hyphy busted --alignment ${file}_uniq.fas --tree ${file}_tree.txt" >> run_busted_${file}.sh
done < "uniq"

# create BUSTED batches
ls *busted*.sh > busted_scripts
split --numeric-suffixes=1 -l 4000 busted_scripts busted_batch

while IFS= read -r file
do
  sbatch $file
done < "busted_batch01"

while IFS= read -r file
do
  sbatch $file
done < "busted_batch02"

while IFS= read -r file
do
  sbatch $file
done < "busted_batch03"

## make aBSREL job scripts

# make scripts for runs without duplicate seqs
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> run_absrel_${file}.sh
  echo -e "#SBATCH -o logs/${file}_absrel.out" >> run_absrel_${file}.sh
  echo -e "#SBATCH -e logs/${file}_absrel.err" >> run_absrel_${file}.sh
  echo -e "#SBATCH -p shared" >> run_absrel_${file}.sh
  echo -e "#SBATCH -n 1" >> run_absrel_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> run_absrel_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> run_absrel_${file}.sh
  echo -e "source activate align\n" >> run_absrel_${file}.sh
  echo -e "hyphy absrel --alignment ${file}_nuc.fa_codon.msa --tree ${file}_tree.txt" >> run_absrel_${file}.sh
done < "split_aligns"

# make scripts for runs with duplicate seqs 
while IFS= read -r file
do
  echo -e '#!/bin/bash' >> run_absrel_${file}.sh
  echo -e "#SBATCH -o logs/${file}_absrel.out" >> run_absrel_${file}.sh
  echo -e "#SBATCH -e logs/${file}_absrel.err" >> run_absrel_${file}.sh
  echo -e "#SBATCH -p shared" >> run_absrel_${file}.sh
  echo -e "#SBATCH -n 1" >> run_absrel_${file}.sh
  echo -e "#SBATCH -t 48:00:00" >> run_absrel_${file}.sh
  echo -e "#SBATCH --mem=9000\n" >> run_absrel_${file}.sh
  echo -e "source activate align\n" >> run_absrel_${file}.sh
  echo -e "hyphy absrel --alignment ${file}_uniq.fas --tree ${file}_tree.txt" >> run_absrel_${file}.sh
done < "uniq"

# create aBSREL batches
ls *absrel*.sh > absrel_scripts
split --numeric-suffixes=1 -l 4000 absrel_scripts absrel_batch

while IFS= read -r file
do
  sbatch $file
done < "absrel_batch01"

while IFS= read -r file
do
  sbatch $file
done < "absrel_batch02"

while IFS= read -r file
do
  sbatch $file
done < "absrel_batch03"
