#################
## OrthoFinder ##
#################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos

#conda create -c bioconda -n ortho orthofinder gffread
conda activate ortho


# copy GTFs from CESAR 
mkdir -p from_cesar/input_data
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_* input_data/
cd from_cesar/
# translate to protein seq
# NB: this step creates an an index file for the reference genome in the 2bitdir - if you don't want that to happen, copy the fastas, then run in the input_data dir
gffread cleaned_reordered_hetAtr.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/hetAtr/hetAtr.fasta -y hetAtr.translated.fa -S 
gffread cleaned_reordered_netAur.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/netAur/netAur.fasta -y netAur.translated.fa -S
gffread cleaned_reordered_oxyJam.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/oxyJam/oxyJam.fasta -y oxyJam.translated.fa -S
gffread cleaned_reordered_stiNae.sorted.gtf -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/stiNae/stiNae.fasta -y stiNae.translated.fa -S
# tidy
mv *.translated.fa input_data/

# get the rest of the species GFFs from the WGA 
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

# convert the Genbank formats to GFF
for file in *.gbff;
do
  mv $file $file.gb
  /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/brename -p ".gbff.gb" -r ".gb"
  singularity exec /n/singularity_images/informatics/maker/maker\:2.31.10--pl526_16.sif bp_genbank2gff.pl --stdout --file $file --source genbank > $file.gff
  /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/brename -p ".gb.gff" -r ".gff"
done

# rename the NCBI GFFs 
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff ansCyg.gff
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff cotJap.gff
mv GCF_002078875.1_NumMel1.0_genomic.gff numMel.gff
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff anaPla.gff

# translate to protein seq
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal numMel syrMik tymCupPin;
do
  gffread $file.gff -g /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/$file/$file.fasta -y $file.translated.fa -S
done

# had indexing issues with galGal, need to strip version information from GFF
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
