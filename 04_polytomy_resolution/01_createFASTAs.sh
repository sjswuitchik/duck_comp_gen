# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# copy over CNEE data
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_cnees/hetAtr/galGal6_final_merged_CNEEs_named.bed  /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_cnees/hetAtr/all_cnees.tab .

# create keys from from full list of CNEEs
#conda create -n r -c bioconda r-base r-tidyverse
source activate r
Rscript cnee_keys.R

# write one FASTA per CNEE
awk 'NR == FNR { cnee[$1] } NR != FNR && $2 in cnee { printf(">%s\n%s\n", $1, $3) >> $2 ".fa"; close($2 ".fa"); }' spp_by_cnee_allSpp.tsv  all_cnees.tab
