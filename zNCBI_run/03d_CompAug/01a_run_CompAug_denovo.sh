## in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes/

module load Anaconda/5.0.1-fasrc01 samtools/1.10-fasrc01
#conda create -c conda-forge -c bioconda -n compAug augustus star
source activate compAug

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

# rename 
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.fna.gz anaPla.ncbi.fasta.gz
mv GCA_002592135.1_ASM259213v1_genomic.fna.gz ansBra.ncbi.fasta.gz
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna.gz ansCyg.ncbi.fasta.gz
mv GCA_006229135.1_ASM622913v1_genomic.fna.gz ansInd.ncbi.fasta.gz
mv GCA_006130075.1_GSC_cangoose_1.0_genomic.fna.gz braCan.ncbi.fasta.gz
mv GCA_008692595.1_Cv_LA_1.0_genomic.fna.gz colVir.ncbi.fasta.gz
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.fna.gz gotJap.ncbi.fasta.gz
mv GCF_000002315.6_GRCg6a_genomic.fna.gz galGal.ncbi.fasta.gz
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz hetAtr.ncbi.fasta.gz
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz netAur.ncbi.fasta.gz
mv GCF_002078875.1_NumMel1.0_genomic.fna.gz numMel.ncbi.fasta.gz
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.ncbi.fasta.gz
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz stiNae.ncbi.fasta.gz
mv GCA_003435085.1_NTU_Smik_1.2_genomic.fna.gz syrMik.ncbi.fasta.gz
mv GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1_genomic.fna.gz tymCupPin.ncbi.fasta.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/075/105/GCA_011075105.1_BPBGC_Hatr_1.0/GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/074/415/GCA_011074415.1_BPBGC_Snae_1.0/GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/077/185/GCA_011077185.1_BPBGC_Ojam_1.0/GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/076/525/GCA_011076525.1_BPBGC_Naur_1.0/GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz
 
mv GCA_011075105.1_BPBGC_Hatr_1.0_genomic.fna.gz hetAtr.ncbi.fasta.gz
mv GCA_011074415.1_BPBGC_Snae_1.0_genomic.fna.gz stiNae.ncbi.fasta.gz
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.ncbi.fasta.gz
mv GCA_011076525.1_BPBGC_Naur_1.0_genomic.fna.gz netAur.ncbi.fasta.gz

gunzip *.gz
cd ..

# load genomes into an SQLite database
while read line
do
  species=$(echo "$line" | cut -f 1)
  genome=$(echo "$line" | cut -f 2)
  load2sqlitedb --noIdx --species=$species --dbaccess=chicken.db $genome
done <genomes.tbl

# index
load2sqlitedb --makeIdx --dbaccess=chicken.db

# qc
sqlite3 -header -column chicken.db "\
 SELECT speciesname, \
  sum(end-start+1) AS 'genome length',\
  count(*) AS '# chunks',\
  count(distinct seqnr) AS '# seqs'\
 FROM genomes natural join speciesnames\
 GROUP BY speciesname;"

# set up de novo run - first convert HAL to MAF with comp Aug specific conversion
mkdir augCGP_denovo
cd augCGP_denovo
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03b_cesar/gallo_ncbi.hal .

# separate HAL into chunks for more efficient processing
mkdir mafs
# chunk size and overlap are recommendations for WGAs
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif #NB could also use the augustus image augustus-2020-05-27-1b69b25ed001.sif
hal2maf_split.pl --halfile gallo_ncbi.hal --refGenome galGal --cpus 8 --chunksize 2500000 --overlap 500000 --outdir mafs

# assign numbers to alignment chunks
num=1
for f in mafs/*.maf; 
do 
 ln -s $f $num.maf; ((num++)); 
done

# run Comp Aug de novo without hints
for ali in mafs/*.maf;
do
id=${ali%.maf} 
augustus \
--species=chicken \
--softmasking=1 \
--treefile=top1.nwk \
--alnfile=$ali \
--dbaccess=../chicken.db \
--speciesfilenames=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/genomes.tbl \
--alternatives-from-evidence=0 \
--/CompPred/outdir=pred$id > aug$id.out 2> err$id.out &
done
