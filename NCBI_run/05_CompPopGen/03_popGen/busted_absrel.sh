# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy emboss
source activate busted

# bring over protein FASTAs and SeqIDs from OrthoFinder
mkdir fastas
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/WorkingDirectory/Species*.fa /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/WorkingDirectory/SequenceIDs.txt fastas/
for file in fastas/*.fa;
do
  awk 'NR == FNR {seqid[">"substr($1, 1, length($1)-1)] = $2; next} /^>/ { print ">" seqid[$1]; next} {print}' fastas/SequenceIDs.txt $file > $file.named
done



