# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy emboss
source activate busted

mkdir single_orthos
mkdir gene_trees
mkdir translated
# bring over protein alignments and gene trees from OrthoFinder
# arg list too long to do a regular cp & working on holylfs is a nightmare, so do this instead to get things over to holyscratch
cd /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/Single_Copy_Orthologue_Sequences
for file in *.fa;
do
  cp -v $file /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/single_orthos/
done

cd ../
cp Orthogroups/Orthogroups_SingleCopyOrthologues.txt Gene_Trees/
cd Gene_Trees

sed 's/$/_tree.txt/' Orthogroups_SingleCopyOrthologues.txt > tmp && mv tmp Orthogroups_SingleCopyOrthologues.txt

for file in `cat Orthogroups_SingleCopyOrthologues.txt`;
do
  cp -v $file /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/gene_trees
done

# translate gene names to species names
cd /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted/
./translate.sh

# backtranslate protein seq to DNA
cd translate/
for file in *.fa;
do
  backtranseq -sequence $file -sprotein1 -auto -cfile Echick.cut -outfile $file.dna
done

# align with PRANK
sbatch run_prank.sh

# run HmmCleaner on PRANK alignments
for file in *.fas;
do
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl $file
done



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
