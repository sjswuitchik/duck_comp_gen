#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos

#conda create -c bioconda -n ortho orthofinder gffread
conda activate ortho

# copy GTFs from CESAR 
mkdir -p from_cesar/input_data
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_* from_cesar/
cd from_cesar/
# translate to protein seq
# NB: this step creates an an index file for the reference genome in the 2bitdir - if you don't want that, copy the fastas over, then run in your working dir
for file in hetAtr netAur oxyJam stiNae;
do
  gffread cleaned_reordered_$file.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta -y $file.translated.fa -S 
  mv $file.translated.fa input_data/
done

# get the rest of the species annotations from the WGA 
cd ../
mkdir -p from_ncbi/input_data
cd from_ncbi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/850/225/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/971/095/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/078/875/GCF_002078875.1_NumMel1.0/GCF_002078875.1_NumMel1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/592/135/GCA_002592135.1_ASM259213v1/GCA_002592135.1_ASM259213v1_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/229/135/GCA_006229135.1_ASM622913v1/GCA_006229135.1_ASM622913v1_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/130/075/GCA_006130075.1_GSC_cangoose_1.0/GCA_006130075.1_GSC_cangoose_1.0_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/692/595/GCA_008692595.1_Cv_LA_1.0/GCA_008692595.1_Cv_LA_1.0_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/435/085/GCA_003435085.1_NTU_Smik_1.2/GCA_003435085.1_NTU_Smik_1.2_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/870/855/GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1/GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1_genomic.gbff.gz

gunzip *.gz

mv GCA_002592135.1_ASM259213v1_genomic.gbff ansBra.gb
mv GCA_006229135.1_ASM622913v1_genomic.gbff ansInd.gb
mv GCA_006130075.1_GSC_cangoose_1.0_genomic.gbff braCan.gb
mv GCA_008692595.1_Cv_LA_1.0_genomic.gbff colVir.gb
mv GCA_003435085.1_NTU_Smik_1.2_genomic.gbff syrMik.gb
mv GCA_001870855.1_T_cupido_pinnatus_GPC_3440_v1_genomic.gbff tymCupPin.gb

# convert the Genbank formats to GFF
for file in ansBra ansInd braCan colVir syrMik tymCupPin;
do
  singularity exec /n/singularity_images/informatics/maker/maker\:2.31.10--pl526_16.sif bp_genbank2gff.pl --stdout --file $file.gb --source genbank > $file.gff
done

# rename the NCBI GFFs 
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff ansCyg.gff
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff cotJap.gff
mv GCF_002078875.1_NumMel1.0_genomic.gff numMel.gff
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff anaPla.gff
mv GCF_000002315.6_GRCg6a_genomic.gff galGal.gff

# translate to protein seq
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap numMel syrMik tymCupPin;
do
  gffread $file.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta -y $file.translated.fa -S
done

# need to strip version information from galGal GFF
python3 stripGFF.py galGal.gff > galGal.stripped.gff
gffread galGal.stripped.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/galGal/galGal.fasta -y galGal.translated.fa -S

# and conversion issues with ansBra, ansInd, braCan, colVir, syrMik, tymCupPin; need to strip version information from FASTA
for file in ansBra ansInd braCan colVir syrMik tymCupPin;
do
  cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta .
  awk -f stripFasta.awk $file.fasta > $file.stripped.fasta
  gffread $file.gff -g $file.stripped.fasta -y $file.translated.fa -S 
done

# tidy
mv *.translated.fa input_data/

# run OrthoFinder 
cd ..
mkdir -p run_ortho/input_data
cp -v from_cesar/input_data/* from_ncbi/input_data/* run_ortho/input_data
orthofinder -f run_ortho/input_data
