# currently a janky mix of manipulation in R and awk
conda activate r

Rscript _____.R 

awk 'BEGIN { FS = OFS = ","} NR == 1 {print "file,total_branches,sig_branches,hetAtr_gene"; next} { match($0, /hetAtr_[^,]+/); hetAtr_gene = substr($0, RSTART, RLENGTH); print $1, $2, $3, hetAtr_gene}' abs_cleanR.csv > abs_clean_final.csv

Rscript ____.R
# download Entrez ID records from https://www.ncbi.nlm.nih.gov/sites/batchentrez using galGal_protIDs.tsv, and 'Send To' File > GenPept with 'Show GI' option
./parse_gp.awk sequence.gp > entrezIDs.tsv

# download the background records from NCBI Protein with same 'Send to' options
./parse_gp.awk sequence.gp > bg_entrezIDs.tsv

Rscript ____.R
