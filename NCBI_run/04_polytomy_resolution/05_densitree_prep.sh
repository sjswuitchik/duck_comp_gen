# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/trimmed/finished_fastas/

mkdir densitree_input
for file in *.bestTree;
do
  cp $file densitree_input/
done

cd densitree_input/

for file in *.bestTree;
do
  python3 newick2nexus.py -i $file -o $file.nex
done

