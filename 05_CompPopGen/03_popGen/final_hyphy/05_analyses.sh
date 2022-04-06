
conda activate r

Rscript _____.R 

# download Entrez ID records from https://www.ncbi.nlm.nih.gov/sites/batchentrez using galGal_protIDs.tsv, and 'Send To' File > GenPept with 'Show GI' option
./parse_gp.awk sequence.gp > entrezIDs.tsv

# download the background records from NCBI Protein with same 'Send to' options
./parse_gp.awk sequence.gp > bg_entrezIDs.tsv

Rscript ____.R
