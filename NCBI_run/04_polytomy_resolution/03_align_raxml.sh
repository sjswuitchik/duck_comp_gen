# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# submit both MAFFT and MUSCLE alignments 
sbatch mafftAlign.sh
sbatch muscleAlign.sh

##### CNEE alignments from MUSCLE were manually trimmed to remove sites where indels where present in >1 spp

# rename CNEE alignments to deal with that one stupid parsing mistake from way back
#conda create -n raxml -c bioconda raxml-ng modeltest-ng rename 
source activate raxml 
for file in trimmed/*.fa;
do
  rename -d zf $file
done

# model selection for RAxML-NG
mkdir -p trimmed/subset
cd trimmed/
for file in $(ls -p | grep -v / | tail -50)
do
  cp $file subset/
done

for file in subset/*.fa;
do
  modeltest-ng -d nt -i $file -o $file -T raxml
done

for file in subset/*.fa;
do
  raxml-ng --check --msa $file --model HKY+G4 --prefix tcheck_$file
done

cd ../
sbatch run_raxml.sh

# test --all against simple inference in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/trimmed/finished_fastas
sbatch run_raxml_all.sh


# protein coding alignments in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/raxml_prot_code
# copy over alignments from /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft/clean_align/all_spp & /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft/last2hmmClean
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft/clean_align/all_spp/*.filtered .
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/aligned/mafft/last2hmmClean/*.filtered .
