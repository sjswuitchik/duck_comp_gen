# in /n/holyscratch01/informatics/swuitchik/ducks/snakemake/hetAtr_stiNae_qc

source activate vcfqc

# download PRANK
wget http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
tar zxvf http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
rm prank.linux64.170427.tgz 
# grab gallo HAL from WGA
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/02_wga/gallo_ncbi.hal .
