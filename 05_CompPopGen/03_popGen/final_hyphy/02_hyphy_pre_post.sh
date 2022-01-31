## in /n/holyscratch01/informatics/swuitchik/ducks/compGen/

conda activate align

# copy over gene trees
mkdir gene_trees
cut -f1 /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/clean_orthos.tsv | sed -e "1d" > clean_ogs

while IFS= read -r file
do
  cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Gene_Trees/${file}_tree.txt gene_trees/
done < "clean_ogs"

# prep hyphy and align
cd og_fastas/
git clone https://github.com/veg/hyphy-analyses.git
sbatch run_hyphy_prep.sh

# check for missed alignments
ls *_codon.msa | sort > codon_aligns
ls *_nuc.fa | sort > og_aligns
sed -i 's/\_codon\.msa//g' codon_aligns
sed -i 's/\_nuc\.fa//g' og_aligns
comm -3 codon_aligns og_aligns 
#OG0003606
#OG0008416
#OG0009185
#OG0010355
#OG0013409
#OG0013485

# redo run_hyphy_prep.sh for missed groups - screen 10549

# clean alignments
sbatch run_hmmCleaner.sh

# reorg dir
mkdir inter_files
mv *_copies.json *_filtered.json *_nuc.fa_nuc.fas *_nuc.fa_protein.fas *_protein.msa *_hmm.score *_hmm.log inter_files/

# fix gene trees to match alignment outputs
cp -vr ../gene_trees/*.txt .

for file in *.txt;
do
  sed -i 's/_\([[:alnum:]]*\)\./_\1_/g' $file
#  cp -v $file hyphy-analyses/remove-duplicates/
done

# remove duplicate sequences and trim tree accordingly, outputting FASTA rather than NH
#while IFS= read -r file
#do
#  hyphy hyphy-analyses/remove-duplicates/remove-duplicates.bf --msa ${file}_nuc.fa_codon.msa --tree ${file}_tree.txt --output ${file}_uniq.fas ENV="DATA_FILE_PRINT_FORMAT=9"
#done < "../clean_ogs.tsv"

# check for uniq vs dup
#ls *_uniq.fas > uniq
#sed -i 's/\_uniq\.fas//g' uniq
#comm -3 codon_aligns uniq > split_aligns

ls *_hmm.fasta | sort > clean_aligns
sed -i 's/\_codon\_hmm\.fasta//g' clean_aligns
comm -3 codon_aligns clean_aligns
# missing OG0003786 from clean_aligns, re-run hmmCleaner on it, rm clean_aligns and re-make
ls *_hmm.fasta | sort > clean_aligns
sed -i 's/\_codon\_hmm\.fasta//g' clean_aligns

# prune trees 
while IFS= read -r file
do
  grep '^>' ${file}_codon.msa > ${file}_tips
  sed -i 's/>//g' ${file}_tips 
  nw_prune -v -f ${file}_tree.txt ${file}_tips > ${file}_prunedTree.txt
done < "clean_aligns"
