# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# create env with both alignment programs
#conda create -n align -c bioconda mafft muscle
conda activate align

# run subset of 50 FASTAs through MUSCLE & MAFFT
mkdir subset
cd fastas/
for file in $(ls -p | grep -v / | tail -50)
do
  mv $file ../subset/
done
cd ..

for file in subset/*.fa;
do
  muscle -in $file -quiet -out $file.afa
done 

for file in subset/*.fa;
do
  mafft --quiet $file > $file.mafft
done
