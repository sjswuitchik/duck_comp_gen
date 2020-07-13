
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/longranger_proc

module load longranger/2.2.2-fasrc01

# make references from fastas
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/077/185/GCF_011077185.1_BPBGC_Ojam_1.0/GCF_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
gunzip *.gz

longranger mkref GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna
longranger mkref GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna
longranger mkref GCF_011077185.1_BPBGC_Ojam_1.0_genomic.fna
longranger mkref GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna

# run longranger
longranger wgs --id hetAtr \
--fastqs=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/longranger_proc/raw_fastqs/04_hetAtr \
--sex=female \
--reference=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/longranger_proc/refdata-GCA_011075105.1_BPBGC_Hatr_1.0_genomic \
--vcmode=gatk:/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar

