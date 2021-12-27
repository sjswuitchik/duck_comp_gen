## in /n/holyscratch01/informatics/swuitchik/ducks/compGen/

git clone https://github.com/veg/hyphy-analyses.git
chmod +x hyphy-analyses/remove-duplicates/remove-duplicates.bf 

# copy over gene trees
mkdir gene_trees
cut -f1 /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/clean_orthos.tsv | sed -e "1d" > clean_ogs.tsv

while IFS= read -r file
do
  cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Gene_Trees/${file}_tree.txt gene_trees/
done < "clean_ogs.tsv"

conda activate align

# prep hyphy and align
cd og_fastas/
sbatch run_hyphy_prep.sh

# check for missed alignments
ls *_codon.msa | sort > uniq_aligns
ls *_nuc.fa | sort > og_aligns
sed -i 's/\_nuc\.fa\_codon\.msa//g' uniq_aligns
sed -i 's/\_nuc\.fa//g' og_aligns
comm -3 uniq_aligns og_aligns 
#OG0003606
#OG0008416
#OG0009185
#OG0010355
#OG0013409
#OG0013485
# redo run_hyphy_prep.sh for missed groups

# fix gene trees to match alignment outputs
cp -vr ../gene_trees/*.txt .

for file in *.txt;
do
  sed -i 's/_\([[:alnum:]]*\)\./_\1_/g' $file
  cp -v $file hyphy-analyses/remove-duplicates/
done

# remove duplicate sequences and trim tree accordingly
while IFS= read -r file
do
  hyphy hyphy-analyses/remove-duplicates/remove-duplicates.bf --msa ${file}_nuc.fa_codon.msa --tree ${file}_tree.txt --output ${file}_uniq.nh
done < "../clean_ogs.tsv"

# check for uniq vs dup
ls *_uniq.nh > uniq
sed -i 's/\_uniq\.nh//g' uniq
comm -3 uniq_aligns uniq > split_aligns

