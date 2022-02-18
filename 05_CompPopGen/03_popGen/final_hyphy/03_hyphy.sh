# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/og_fastas

# make BUSTED job scripts
mkdir logs

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
  echo -e "hyphy busted --alignment ${file}_codon_hmm.fasta --tree ${file}_prunedTree.txt" >> run_busted_${file}.sh
done < "clean_aligns"

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

# make aBSREL job scripts
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
  echo -e "hyphy absrel --alignment ${file}_codon_hmm.fasta --tree ${file}_prunedTree.txt" >> run_absrel_${file}.sh
done < "clean_aligns"

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

## missing some BUSTED outputs, so resubmit 
ls -lh *BUSTED.json | cut -c 24- | sort -nr | awk '$2 == 0 {print $6}' | sort > failed_busted
sed -i 's/\_codon\_hmm\.fasta\.BUSTED\.json//g' failed_busted 
ls *BUSTED.json > all_bust
sed -i 's/\_codon\_hmm\.fasta\.BUSTED\.json//g' all_bust
comm -3 clean_aligns all_bust > rerun_busted 

while IFS= read -r file
do
  sbatch run_busted_${file}.sh
done < "rerun_busted"

## some aBSREL submissions were interrupted, resub
ls *ABSREL.json > all_abs
sed -i 's/\_codon\_hmm\.fasta\.ABSREL\.json//g' all_abs
comm -3 clean_aligns all_abs > rerun_abs
comm -3 rerun_busted rerun_abs
# same alignments failed both, looks like it's a dup seq issue
# remake scripts with new workflow for dup seqs 

while IFS= read -r file
do
  rm run_busted_${file}.sh
  rm run_absrel_${file}.sh
done < "rerun_busted"

while IFS= read -r file
do
  hyphy hyphy-analyses/remove-duplicates/remove-duplicates.bf --msa ${file}_codon.msa --output ${file}_uniq.fa ENV="DATA_FILE_PRINT_FORMAT=9"
done < "rerun_busted"

while IFS= read -r file
do
  hyphy hyphy-analyses/remove-duplicates/remove-duplicates.bf --msa ${file}_codon.msa --output ${file}_uniq.fa ENV="DATA_FILE_PRINT_FORMAT=9"
done < "rerun_abs"

ls *_uniq.fa > removed_dups
sed -i 's/\_uniq\.fa//g' removed_dups
comm -3 rerun_busted removed_dups 
# every entry the same

while IFS= read -r file
do
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl ${file}_uniq.fa
  grep '^>' ${file}_uniq_hmm.fasta > ${file}_tips
  sed -i 's/>//g' ${file}_tips
  nw_prune -v -f ${file}_tree.txt ${file}_tips > ${file}_prunedTree.txt
done < "removed_dups"

while IFS= read -r file
do
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl ${file}_uniq.fa
  grep '^>' ${file}_uniq_hmm.fasta > ${file}_tips
  sed -i 's/>//g' ${file}_tips
  nw_prune -v -f ${file}_tree.txt ${file}_tips > ${file}_prunedTree.txt
done < "rerun_abs"

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
  echo -e "hyphy busted --alignment ${file}_uniq_hmm.fasta --tree ${file}_prunedTree.txt" >> run_busted_${file}.sh
done < "rerun_busted"

while IFS= read -r file
do
  sbatch run_busted_${file}.sh
  sleep 2
done < "rerun_busted"

ls *BUSTED.json > all_bus
sed -i 's/\_codon\_hmm\.fasta\.BUSTED\.json//g' all_bus
sed -i 's/\_uniq\_hmm\.fasta\.BUSTED\.json//g' all_bus
comm -3 clean_aligns all_bus > rerun_bus

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
  echo -e "hyphy absrel --alignment ${file}_uniq_hmm.fasta --tree ${file}_prunedTree.txt" >> run_absrel_${file}.sh
done < "rerun_abs"

while IFS= read -r file
do
  sbatch run_absrel_${file}.sh
  sleep 2
done < "rerun_abs"
