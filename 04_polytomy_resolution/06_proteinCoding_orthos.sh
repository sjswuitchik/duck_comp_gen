# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_coding
module load jdk/16-fasrc01
mkdir orthos

# copy over orthogroup alignments from HyPhy runs
cp -v ../compGen/og_fastas/*_codon_hmm.fasta ../compGen/og_fastas/*_uniq_hmm.fasta orthos/

# run RAxML
sbatch run_raxml_prot.sh

# once RAxML has finished
cat orthos/*.bestTree > prot.tree

cat orthos/*.support > prot.support.tree

sed 's/_[^:]*//g' prot.tree > prot.trim.tee

git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
./make.sh
cp ../prot.trim.tree .

# run ASTRAL
sbatch run_astral_prot.sh
