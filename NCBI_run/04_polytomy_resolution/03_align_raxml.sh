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

# combine best trees into ASTRAL input file
cd subset/
cat *.bestTree > test.tree
