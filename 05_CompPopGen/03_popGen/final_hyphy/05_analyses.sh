# currently a janky mix of manipulation in R and awk
conda activate r

Rscript clean_absrel.R 

awk 'BEGIN { FS = OFS = ","} NR == 1 {print "file,total_branches,sig_branches,hetAtr_gene"; next} { match($0, /hetAtr_[^,]+/); hetAtr_gene = substr($0, RSTART, RLENGTH); print $1, $2, $3, hetAtr_gene}' abs_cleanR.csv > abs_clean_final.csv

Rscript absrel_ortho_wrangling.R

# download Entrez ID records from https://www.ncbi.nlm.nih.gov/sites/batchentrez using galGal_protIDs.tsv, and 'Send To' File > GenPept with 'Show GI' option
./parse_gp.awk sequence.gp > entrezIDs.tsv

split -l 3000 all_galGal_protIDs.tsv all_galGal_prot

for file in *.gp;
do
  ./parse_gp.awk $file > ${file}.tsv
done

# download the background records from NCBI datasets (gbff) 
./parse_gp.awk genomic.gbff > bg_entrezIDs.tsv

Rscript go_analyses.R
