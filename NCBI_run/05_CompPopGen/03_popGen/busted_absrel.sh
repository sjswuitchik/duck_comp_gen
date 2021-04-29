# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy emboss
source activate busted

# bring over protein alignments from OrthoFinder
mkdir align
# arg too long, work in OrthoFinder working dir for initial translation stuff - come back to this after meeting
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/Orthogroup_Sequences/*.fa align/

# translate the def lines from OrthoFinder to actual seq IDs
for file in fastas/*.fa;
do
  awk 'NR == FNR {seqid[">"substr($1, 1, length($1)-1)] = $2; next} /^>/ { print ">" seqid[$1]; next} {print}' fastas/SequenceIDs.txt $file > $file.named
done

# rename the species FASTAs to the actual species codes
while read id file; 
do
  mv Species${id%:}.fa.named ${file%%.translated.fa}.fa;
done < SpeciesIDs.txt

# tidy
mkdir orig_fastas
mv Species*.fa orig_fastas

