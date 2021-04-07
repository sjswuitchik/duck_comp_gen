# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee
module load jdk/10.0.1-fasrc01

# download ASTRAL
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.7.zip
unzip Astral.5.7.7.zip
rm Astral.5.7.7.zip
cd Astral

# copy data over 
cp ../trimmed/subset/final.tree

# run ASTRAL 
java -jar astral.5.7.7.jar -i final.tree -o final.astral.tree
