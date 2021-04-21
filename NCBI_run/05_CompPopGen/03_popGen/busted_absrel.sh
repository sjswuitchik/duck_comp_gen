# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

source activate vcfqc

# download PRANK
wget http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
tar zxvf http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
rm prank.linux64.170427.tgz 

# grab gallo HAL from WGA
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/02_wga/gallo_ncbi.hal .
# convert to MAF
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif
hal2maf gallo_ncbi.hal gallo_ncbi.maf --refGenome galGal --noAncestors --noDupes
mafTools/bin/mafValidator.py --maf=gallo_ncbi.maf


