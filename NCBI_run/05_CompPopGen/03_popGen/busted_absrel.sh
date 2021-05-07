# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy emboss perl perl-app-cpanminus
source activate busted

# bring over protein alignments from OrthoFinder
mkdir align
# arg list too long to do a regular cp & working on holylfs is a nightmare, so do this instead to get things over to holyscratch
cd /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/Single_Copy_Orthologue_Sequences
for file in *.fa;
do
  cp -v $file /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/align/
done
cd /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/align/

# backtranslate protein seq to DNA
for file in *.fa;
do
  backtranseq -sequence $file -sprotein1 -auto -cfile Echick.cut -outfile $file.dna
done

# align with PRANK
sbatch run_prank.sh

# need to containerize HmmCleaner - nasty list of perl dependencies

# run HmmCleaner on PRANK alignments
ls *.fas > inList
HmmCleaner.pl inList 2> hmm.err


############# 

## translate the def lines from OrthoFinder to actual seq IDs
#for file in fastas/*.fa;
#do
#  awk 'NR == FNR {seqid[">"substr($1, 1, length($1)-1)] = $2; next} /^>/ { print ">" seqid[$1]; next} {print}' fastas/SequenceIDs.txt $file > $file.named
#done
#
## rename the species FASTAs to the actual species codes
#while read id file; 
#do
#  mv Species${id%:}.fa.named ${file%%.translated.fa}.fa;
#done < SpeciesIDs.txt
#
## tidy
#mkdir orig_fastas
#mv Species*.fa orig_fastas
