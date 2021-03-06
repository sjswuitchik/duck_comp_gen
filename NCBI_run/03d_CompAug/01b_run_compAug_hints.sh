#### run Comp Aug with RNA seq hints from ruddy duck (oxyJam) and hints from NCBI GFFs of chicken (galGal6), mallard (anaPla), and ruddy duck (oxyJam) ####
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug
## from http://bioinf.uni-greifswald.de/augustus/binaries/tutorial-cgp/rnaseq.html and http://bioinf.uni-greifswald.de/augustus/binaries/tutorial-cgp/cactus.html#runCactus

module load Anaconda/5.0.1-fasrc01 samtools/1.10-fasrc01 perl/5.26.1-fasrc01
#conda create -c conda-forge -c bioconda -n compAug augustus star
source activate compAug

mkdir -p 03d_CompAug/augCGP_rnahints/
cd augCGP_rnahints/

# get oxyJam genome & annotations for oxyJam, galGal, anaPla
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/077/185/GCF_011077185.1_BPBGC_Ojam_1.0/GCF_011077185.1_BPBGC_Ojam_1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/077/185/GCA_011077185.1_BPBGC_Ojam_1.0/GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/850/225/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5/GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff.gz
mv GCF_011077185.1_BPBGC_Ojam_1.0_genomic.gff.gz oxyJam.gff.gz
mv GCA_011077185.1_BPBGC_Ojam_1.0_genomic.fna.gz oxyJam.fa.gz
mv GCF_000002315.6_GRCg6a_genomic.gff.gz galGal.gff.gz
mv GCF_003850225.1_IASCAAS_PekingDuck_PBH1.5_genomic.gff.gz anaPla.gff.gz

mkdir genome
mv oxyJam.fa genome/
mkdir gffs
mv oxyJam.gff.gz galGal.gff.gz anaPla.gff.gz gffs/
gunzip genome/* gffs/*

# get RNA seq fastqs from SRA (Ruddy Duck replicates 1-3 from bioProject PRJNA517454) 
mkdir fastqs/
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz
tar zxvf sratoolkit.2.10.8-centos_linux64.tar.gz
rm *.tar.gz
mv sratoolkit.2.10.8-centos_linux64/ sratoolkit/
cd sratoolkit/bin/
./fasterq-dump -O ../../fastqs/ -p -S --skip-technical SRR8495655 SRR8495656 SRR8495657 SRR8495658 SRR8495659 SRR8495660 SRR8495661 SRR8495662 SRR8495663 SRR8495664 SRR8495665 SRR8495666 SRR8495667 SRR8495668 SRR8495669 SRR8495670 SRR8495641 SRR8495642 SRR8495643 SRR8495644 SRR8495645 SRR8495646 SRR8495647 SRR8495648 SRR8495649 SRR8495650 SRR8495651 SRR8495652 SRR8495653 SRR8495654 SRR8495625 SRR8495626 SRR8495627 SRR8495628 SRR8495629 SRR8495630 SRR8495631 SRR8495632 SRR8495633 SRR8495634 SRR8495635 SRR8495636 SRR8495637 SRR8495638 SRR8495639 SRR8495640
cd ../..

# generate genome indices
sbatch run_starPrep.sh

# map
sbatch run_starMap.sh

# 2nd pass
sbatch run_starMap_2pass.sh

# sort output and convert to BAM
sort -k 1,1 oxyJamPass2Aligned.out.sam > oxyJam.input.sort.sam
samtools view oxyJam.input.sort.sam -b -h > oxyJam.input.sort.bam

mkdir star_outputs
mv oxyJamPass* star2* star_7* starprep* star_outputs

# filter & create hints for oxyJam from RNAseq data
mkdir hints/
cd hints/
cp /n/home_rc/afreedman/software/Augustus/augustus-2020-05-27-1b69b25ed001.sif .
singularity shell --cleanenv augustus-2020-05-27-1b69b25ed001.sif
filterBam --in ../oxyJam.input.sort.bam --out oxyJam.filter.bam --uniq --paired --pairwiseAlignment 
bam2hints --in oxyJam.filter.bam --out oxyJam.hints.gff
exit 
samtools sort -o oxyJam.filter.sorted.bam -O bam oxyJam.filter.bam 
singularity shell --cleanenv augustus-2020-05-27-1b69b25ed001.sif
bam2wig oxyJam.filter.sorted.bam > oxyJam.wig
exit
cat oxyJam.wig | ./wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --src=W --type=ep --radius=4.5 > oxyJam.hints.ep.gff

# concat hints
cat oxyJam.hints.gff oxyJam.hints.ep.gff > oxyJam.combohints.gff

# replace chr names in oxyJam NCBI GFF to match database
cd ../gffs/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/077/185/GCF_011077185.1_BPBGC_Ojam_1.0/GCF_011077185.1_BPBGC_Ojam_1.0_assembly_report.txt
sed 's/\r$//g' GCF_011077185.1_BPBGC_Ojam_1.0_assembly_report.txt | grep -v "^#" | cut -f1,5,7 > oxyJam_chr_key
awk '{print $3, $2}' oxyJam_chr_key > acckey
./replace_chrs.pl acckey oxyJam.gff > oxyJam.repl.gff
cd ..

# create hints for galGal, anaPla, and oxyJam using NCBI GFFs
grep -P "\t(CDS|intron)\t" gffs/oxyJam.repl.gff | cut -f1-8 | perl -pe 's/$/\tsource=M/' >hints/oxyJam.ncbihints.gff
grep -P "\t(CDS|intron)\t" gffs/galGal.gff | cut -f1-8 | perl -pe 's/$/\tsource=M/' >hints/galGal.hints.gff
grep -P "\t(CDS|intron)\t" gffs/anaPla.gff | cut -f1-8 | perl -pe 's/$/\tsource=M/' >hints/anaPla.hints.gff

sort -n -k 4,4 hints/oxyJam.ncbihints.gff | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | ./join_mult_hints.pl >temp
mv temp hints/oxyJam.ncbihints.gff
sort -n -k 4,4 hints/galGal.hints.gff | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl >temp
mv temp hints/galGal.hints.gff
sort -n -k 4,4 hints/anaPla.hints.gff | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl >temp
mv temp hints/anaPla.hints.gff

# concat hints
cd  hints/
cat oxyJam.hints.gff oxyJam.hints.ep.gff oxyJam.ncbihints.gff > oxyJam.combohints.gff
cd ..

# clean up intermediate GFFs
mkdir int_files
mv oxyJam.hints.* oxyJam.ncbihints.gff int_files/
cd ..

# load hints into SQL database
for file in hints/*.gff; 
do 
  echo -ne "$(basename $file .hints.gff)\t$file\n"; 
done >hints.tbl
# manually editing this file to correct small formatting issues, will return to this later

# copy database from parent dir (created in 01a_run_CompAug_denovo.sh) 
cp ../chicken.db .
mv chicken.db chicken_rnaseq.db

# load hints into new database
while read line
do
  species=$(echo "$line" | cut -f 1)
  hints=$(echo "$line" | cut -f 2)
  load2sqlitedb --noIdx --species=$species --dbaccess=chicken_rnaseq.db $hints
done <hints.tbl

# finalize db
load2sqlitedb --makeIdx --dbaccess=chicken_rnaseq.db

# qc
sqlite3 -header -column chicken_rnaseq.db "\
  SELECT count(*) AS '#hints',typename,speciesname\
  FROM (hints as H join featuretypes as F on H.type=F.typeid)\
        natural join speciesnames\
  GROUP BY speciesid,typename;"
  
# prep extrinsic config (from https://github.com/nextgenusfs/augustus/blob/master/config/extrinsic/cgp.extrinsic.cfg) with notes from http://bioinf.uni-greifswald.de/augustus/binaries/tutorial-cgp/rnaseq.html

# copy over Newick tree and HAL from de novo dir for simplicity
cp ../augCGP_denovo/top1.nwk .
cp ../augCGP_denovo/gallo_ncbi.hal .

# separate HAL into chunks for more efficient processing
mkdir mafs
# chunk size and overlap are recommendations for WGAs
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif #NB could also use the augustus image augustus-2020-05-27-1b69b25ed001.sif
hal2maf_split.pl --halfile gallo_ncbi.hal --refGenome galGal --cpus 8 --chunksize 2500000 --overlap 500000 --outdir mafs
exit

# assign numbers to alignment chunks
num=1
for f in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/mafs/*.maf; 
do 
 ln -s $f $num.maf; ((num++)); 
done

# run comp Aug with hints
sbatch --array=1-962 run_compAug_singularity.sh

# combine predictions from parallel runs
mkdir joined_pred
while read line
do
  species=$(echo "$line" | cut -f 1)
  find pred* -name "${species}.cgp.gff" >${species}_gtfs.lst;
  joingenes -f ${species}_gtfs.lst -o joined_pred/$species.gff
done < ../genomes.tbl

mkdir slurm_outs/
mv slurm-* slurm_outs/
