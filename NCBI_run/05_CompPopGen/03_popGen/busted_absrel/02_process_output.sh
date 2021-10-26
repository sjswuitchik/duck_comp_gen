# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

mkdir busted_out
mkdir absrel_out
cp busted/aligned/mafft/clean_align/all_spp/*BUSTED* busted_out
cp busted/aligned/mafft/last2hmmClean/*BUSTED* busted_out
cp busted/aligned/mafft/clean_align/all_spp/*ABSREL* absrel_out
cp busted/aligned/mafft/last2hmmClean/*ABSREL* absrel_out

git clone https://github.com/gwct/core.git
python3 core/hyphy-interface/hyphy_to_csv.py -i busted_out -m busted --overwrite -o busted_output
python3 core/hyphy-interface/hyphy_to_csv.py -i absrel_out -m absrel --overwrite -o absrel_output
