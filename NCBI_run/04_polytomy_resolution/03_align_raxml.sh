# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# submit both MAFFT and MUSCLE alignments 
sbatch mafftAlign.sh
sbatch muscleAlign.sh

# rename CNEE alignments to deal with that one stupid parsing mistake from way back
#conda create -n raxml -c bioconda raxml-ng modeltest-ng rename 
source activate raxml 
for file in fastas/muscle/*.afa;
do
  rename -d zf $file
done

for file in fastas/mafft/*.mafft;
do
  rename -d zf $file
done

# model selection for RAxML-NG
mkdir -p fastas/muscle/subset
cd fastas/muscle/
for file in $(ls -p | grep -v / | tail -50)
do
  cp $file subset/
done

for file in subset/*.fa.afa;
do
  modeltest-ng -d nt -i $file -o $file -T raxml
done

for file in subset/*.afa;
do
  raxml-ng --check --msa $file --model HKY+G4 --prefix tcheck_$file
done
