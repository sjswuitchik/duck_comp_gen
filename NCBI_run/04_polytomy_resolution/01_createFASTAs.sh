# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# copy over CNEE data
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_cnees/hetAtr/galGal6_final_merged_CNEEs_named.bed  /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_cnees/hetAtr/all_cnees.tab .


#conda create -n r -c bioconda r-base r-tidyverse
source activate r
Rscript cnee_keys.R
