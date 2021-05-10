# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/trimmed/finished_fastas/

mkdir densitree_input
cp final.all.tree densitree_input/
cd densitree_input/

python3 newick2nexus.py -i final.all.tree -o final.all.tree.nex
