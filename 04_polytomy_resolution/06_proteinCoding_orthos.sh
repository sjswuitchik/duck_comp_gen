# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_coding
module load jdk/10.0.1-fasrc01
mkdir orthos

# copy over orthogroup alignments from HyPhy runs
cp -v ../compGen/og_fastas/*_codon_hmm.fasta ../compGen/og_fastas/*_uniq_hmm.fasta orthos/

# run RAxML
sbatch run_raxml_prot.sh

# once RAxML has finished
cat orthos/*.bestTree > prot.tree

cat orthos/*.support > prot.support.tree

git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
./make.sh
cp ../prot.tree .

# run ASTRAL
java -jar astral.5.7.8.jar -i prot.tree -o prot.astral.tree 2> prot.log
