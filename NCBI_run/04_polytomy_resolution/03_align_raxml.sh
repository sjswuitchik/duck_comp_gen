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
modeltest-ng -d nt -i fastas/muscle/zfCNEE247578.fa.afa -o CNEE275031 -T raxml -v 
