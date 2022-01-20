# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_coding
module load jdk/10.0.1-fasrc01
mkdir orthos

# copy over orthogroup alignments from HyPhy runs
cp -v ../compGen/transfer/* orthos/

# run RAxML
sbatch run_raxml_prot.sh

cd ../Astral
cp ../polytomy_coding/prot.tree .

# run ASTRAL
java -jar astral.5.7.7.jar -i prot.tree -o prot.astral.tree 2> prot.log
