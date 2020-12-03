###############################
## Prep for & running CACTUS ##
###############################

# download all genomes from NCBI for WGA

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/850/225/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/592/135/GCA_002592135.1_ASM259213v1/GCA_002592135.1_ASM259213v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/971/095/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/229/135/GCA_006229135.1_ASM622913v1/GCA_006229135.1_ASM622913v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/130/075/GCA_006130075.1_GSC_cangoose_1.0/GCA_006130075.1_GSC_cangoose_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/692/595/GCA_008692595.1_Cv_LA_1.0/GCA_008692595.1_Cv_LA_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/078/875/GCF_002078875.1_NumMel1.0/GCF_002078875.1_NumMel1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/077/185/GCA_011077185.1_BPBGC_Ojam_1.0/GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/435/085/GCA_003435085.1_NTU_Smik_1.2/GCA_003435085.1_NTU_Smik_1.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/870/855/GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1/GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1_genomic.fna.gz

gunzip *.gz

# rename 
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.fna anaPla.ncbi.fasta
mv GCA_002592135.1_ASM259213v1_genomic.fna ansBra.ncbi.fasta
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna ansCyg.ncbi.fasta
mv GCA_006229135.1_ASM622913v1_genomic.fna ansInd.ncbi.fasta
mv GCA_006130075.1_GSC_cangoose_1.0_genomic.fna braCan.ncbi.fasta
mv GCA_008692595.1_Cv_LA_1.0_genomic.fna colVir.ncbi.fasta
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.fna gotJap.ncbi.fasta
mv GCF_000002315.6_GRCg6a_genomic.fna galGal.ncbi.fasta
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna hetAtr.ncbi.fasta
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna netAur.ncbi.fasta
mv GCF_002078875.1_NumMel1.0_genomic.fna numMel.ncbi.fasta
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna oxyJam.ncbi.fasta
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna stiNae.ncbi.fasta
mv GCA_003435085.1_NTU_Smik_1.2_genomic.fna syrMik.ncbi.fasta
mv GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1_genomic.fna tymCupPin.ncbi.fasta

# need to modify all deflines before running CACTUS 
for file in *.ncbi.fasta; 
do
	sed '/^>/s/ .*//' $file > $file.bre
done

rm *.ncbi.fasta

/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/CESAR2.0/brename -p ".fasta.bre" -r ".fasta" -R

# Run Cactus
sbatch run_cactus.sh 
