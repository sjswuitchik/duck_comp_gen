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
java -jar astral.5.7.7.jar -i final.tree -o final.astral.tree 2> cnee.log

# copy single copy ortho gene trees over from OrthoFinder output to /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee/orthofinder_geneTrees
sed -i.bak 's/_[^:]\{1,\}:/:/g' *.txt
awk '{print}' *_tree.txt > gene_trees.txt
cp gene_trees ../Astral
cd ../Astral
java -jar astral.5.7.7.jar -i gene_trees.txt -o ortho.tree 2> ortho.log

# run ASTRAL on raxml-all run
cp ../trimmed/finished_fastas/final.all.tree .
java -jar astral.5.7.7.jar -i final.all.tree -o final.all.astral.tree 2> all.log
# and on BS tree
cp ../trimmed/finished_fastas/final.support.tree .
java -jar astral.5.7.7.jar -i final.support.tree -o final.support.astral.tree 2> support.log
# run polytomy test
java -jar astral.5.7.7.jar -i final.all.tree -o final.all.poly.tree -t 10 2> all.poly.log
