## in /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees

module load Anaconda/5.0.1-fasrc01
#conda create -c conda-forge -c bioconda -n compAug augustus bcftools htslib samtools
source activate compAug

for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal numMel syrMik tymCupPin;
do
	cp $file.defline.fasta /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/
done 

cd /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/

gzip *.fasta

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/077/185/GCA_011077185.1_BPBGC_Ojam_1.0/GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
 
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz hetAtr.ncbi.fasta.gz
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz stiNae.ncbi.fasta.gz
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.ncbi.fasta.gz
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz netAur.ncbi.fasta.gz



