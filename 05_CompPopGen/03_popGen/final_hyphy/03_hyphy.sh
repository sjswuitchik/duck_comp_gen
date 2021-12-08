# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

# copy over gene trees
mkdir gene_trees
cut -f1 /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/clean_orthos.tsv | sed -e "1d" > clean_ogs.tsv

while IFS= read -r file
do
  cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Gene_Trees/${file}_tree.txt gene_trees/
done < "clean_ogs.tsv"

# edit gene trees to match alignments
cd gene_trees
for file in *.txt;
do
  sed -i 's/\_protein//g' $file 
  sed -i 's/\_translated//g' $file
done

sbatch run_busted.sh
sbatch run_absrel.sh

