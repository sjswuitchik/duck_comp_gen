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
  hyphy remove-duplicates.bf --msa ${file}_codon.msa --tree hyphy-analyses/remove-duplicates/${file}_tree.txt --output ${file}_uniq.nh
done < "../clean_ogs.tsv"
