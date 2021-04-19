# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee
module load jdk/10.0.1-fasrc01

# download ASTRAL
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.7.zip
unzip Astral.5.7.7.zip
rm Astral.5.7.7.zip
cd Astral

# copy data over 
cp ../trimmed/final.tree .

# run ASTRAL 
java -jar astral.5.7.7.jar -i final.tree -o final.astral.tree

# copy single copy ortho gene trees over from OrthoFinder output to /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/orthofinder_geneTrees
sed -i.bak 's/_[^:]\{1,\}:/:/g' *.txt
awk '{print}' *_tree.txt > gene_trees.txt
cp gene_trees ../Astral
cd ../Astral
java -jar astral.5.7.7.jar -i gene_trees.txt -o ortho.tree 2> log

