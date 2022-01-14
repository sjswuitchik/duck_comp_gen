# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

mkdir busted_out
mkdir absrel_out
cp -v og_fastas/*BUSTED.json busted_out
cp -v og_fastas/*ABSREL.json absrel_out

# parse output
git clone https://github.com/gwct/hyphy-interface.git
python3 hyphy-interface/hyphy_to_csv.py -i busted_out -m busted --overwrite -o busted_output.csv
python3 hyphy-interface/hyphy_to_csv.py -i absrel_out -m absrel --overwrite -o absrel_output.csv

# clean up output (keep run meta data but write out file without header for easier analysis)
grep -v '^#' busted_output.csv | sed 's/\_nuc\.fa\_codon\.msa\.BUSTED\.json//g' | sed 's/\_uniq\.fas\.BUSTED\.json//g' > busted_output_clean.csv
grep -v '^#' absrel_output.csv | sed 's/\_nuc\.fa\_codon\.msa\.BUSTED\.json//g' | sed 's/\_uniq\.fas\.BUSTED\.json//g' > absrel_output_clean.csv

# filter aBSREL results for OGs with hetAtr 
grep 'hetAtr' absrel_output_clean.csv > absrel_output_clean_hetAtr.csv
sed -i '1s/^/file\,branches\,num\_ps\_branches\,ps\_pvals\n/' absrel_output_clean_hetAtr.csv
