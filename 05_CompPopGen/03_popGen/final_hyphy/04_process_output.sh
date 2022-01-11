# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

mkdir busted_out
mkdir absrel_out
cp -v og_fastas/*BUSTED.json busted_out
cp -v og_fastas/*ABSREL.json absrel_out

git clone https://github.com/gwct/core.git
python3 core/hyphy-interface/hyphy_to_csv.py -i busted_out -m busted --overwrite -o busted_output.csv
python3 core/hyphy-interface/hyphy_to_csv.py -i absrel_out -m absrel --overwrite -o absrel_output.csv
