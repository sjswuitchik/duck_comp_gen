## in /n/holyscratch01/informatics/swuitchik/ducks/compGen/og_fastas

git clone https://github.com/veg/hyphy-analyses.git
chmod +x hyphy-analyses/remove-duplicates/remove-duplicates.bf 

conda activate align

sbatch run_hyphy_prep.sh

cp -vr ../gene_trees/*.txt .

for file in *.txt;
do
  sed -i 's/_\([[:alnum:]]*\)\./_\1_/g' $file
  cp -v $file hyphy-analyses/remove-duplicates/
done

while IFS= read -r file
do
  hyphy remove-duplicates.bf --msa ${file}_codon.msa --tree hyphy-analyses/remove-duplicates/${file}_tree.txt --output ${file}_uniq.nh
done < "clean_ogs.tsv"
