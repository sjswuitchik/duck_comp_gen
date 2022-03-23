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
grep -v '^#' busted_output.csv | sed 's/\_codon\_hmm\.fasta\.BUSTED\.json//g' | sed 's/\_uniq\_hmm\.fasta\.BUSTED\.json//g' > busted_output_clean.csv
grep -v '^#' absrel_output.csv | sed 's/\_codon\_hmm\.fasta\.ABSREL\.json//g' | sed 's/\_uniq\_hmm\.fasta\ABSREL\.json//g' > absrel_output_clean.csv

# parse aBSREL output to associate FDR-corrected pvalue with node 
./species_pval.awk absrel_output_clean.csv > absrel_output_clean_parse.csv

# grab hetAtr only, just in case 
grep 'hetAtr' absrel_output_clean_parse.csv > absrel_output_clean_parse_hetAtr.csv
sed -n 1p absrel_output_clean_parse.csv > abs_head
cat abs_head absrel_output_clean_parse_hetAtr.csv > absrel_output_clean_parse_hetAtr_head.csv
## nb: this is messy and coded in an airport, come back and fix this up later

# add in Rscript for getting to branch01:branch03

awk 'BEGIN { FS = OFS = ","} NR == 1 {print "file,total_branches,sig_branches,hetAtr_gene"; next} { match($0, /hetAtr_[^,]+/); hetAtr_gene = substr($0, RSTART, RLENGTH); print $1, $2, $3, hetAtr_gene}' abs_cleanR.csv > abs_clean_final.csv
