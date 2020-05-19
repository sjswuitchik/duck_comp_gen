##########################################
## Re-doing ducks with NCBI submissions ##
##########################################

# /n/holyscratch01/informatics/swuitchik/ducks_project/cactus_ncbi 
# FASTAs are already soft masked, so starting from CACTUS

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/077/185/GCA_011077185.1_BPBGC_Ojam_1.0/GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
 
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz hetAtr.ncbi.fasta.gz
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz stiNae.ncbi.fasta.gz
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.ncbi.fasta.gz
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz netAur.ncbi.fasta.gz

gunzip *.gz

cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/*.defline.fasta .
rm *.masked.defline.fasta

for file in *.ncbi.fasta; 
do
	sed '/^>/s/ .*//' $file > $file.bre
done

rm *.ncbi.fasta

/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/brename -p ".fasta.bre" -r ".fasta" -R

# remade galloanserae.txt with new paths for gallo_ncbi.txt
# redid cesar batch script for new files 
