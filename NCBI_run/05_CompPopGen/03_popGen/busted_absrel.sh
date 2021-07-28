# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy bedtools r-base r-tidyverse
source activate busted

mkdir fastas
mkdir gffs

# bring over genomes
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/genomes/* fastas/

# get annotations from NCBI
cd gffs/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/078/875/GCF_002078875.1_NumMel1.0/GCF_002078875.1_NumMel1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/476/345/GCF_015476345.1_ZJU1.0/GCF_015476345.1_ZJU1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/971/095/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0/GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz
gunzip *.gz
# rename
mv GCF_000002315.6_GRCg6a_genomic.gff galGal.gff
mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff ansCyg.gff
mv GCF_001577835.2_Coturnix_japonica_2.1_genomic.gff cotJap.gff
mv GCF_002078875.1_NumMel1.0_genomic.gff numMel.gff
mv GCF_015476345.1_ZJU1.0_genomic.gff anaPla.gff

# bring over annotations from Comp Aug
for file in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/$file.gff .
done

## translate Comp Aug annotations to include gene information based on galGal
# create translation file for galGal transcript to gene
grep -v '#' gffs/galGal.gff | awk '{if ($3 == "mRNA") print $0;}' | python3 galGenes_trans.py > transGene.txt
# quickly reformat transGene file
Rscript reformat.R

# create translation file for Comp Aug spp to galGal transcripts
mkdir trans_files
cd trans_files/
for file in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/02_ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/Orthologues/Orthologues_galGal.translated/galGal.translated__v__$file.translated.tsv .
  sed '1d' galGal.translated__v__$file.translated.tsv | cut -f2,3 > $file_trans.tsv
done

cd ../

# add a gene ID to the GTF with chicken-based genes from the translation files
for sp in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  ./augustus-gtf-transcript-to-gene.sh gffs/$sp.gff transGene_final.tsv trans_files/$sp\_trans.tsv > gffs/$sp\_final.gtf
done

# for NCBI annotations, create CDS-only BEDs
for file in galGal ansCyg cotJap numMel anaPla;
do
  column -s, -t < gffs/$file.gff | awk '$3 == "CDS"' > gffs/$file.cds.gff
  awk -f gff2bed.awk gffs/$file.cds.gff > gffs/$file.cds.bed
  cat gffs/$file.cds.bed | python3 genenames.py > gffs/$file.cds.genes.bed
done

# for Comp Aug annotations, create BEDs with gene names
for file in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  cat gffs/$file\_final.gtf | python3 gffs/genenames_compaug.py | grep -v '^#' | grep 'CDS' | cut -f1,4,5,6 > gffs/$file.cds.genes.bed
done

# extract nucleotide sequences associated with genes
for file in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin galGal ansCyg cotJap numMel anaPla;
do
  bedtools getfasta -fi fastas/$file.ncbi.fasta -bed gffs/$file.cds.genes.bed -name -s > $file.out.fa
done

# fix up with bioawk
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
  bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$file"'", $name, $seq}' $file.out.fa >> all_cds.tab
done

# check 
cut -f1,1 all_cds.tab | sort | uniq -c
#607700 anaPla
#164833 ansBra
#391902 ansCyg
#162112 ansInd
#165052 braCan
#144353 colVir
#566194 cotJap
#669824 galGal
#165003 hetAtr
#165121 netAur
#569903 numMel
#171181 oxyJam
#159217 stiNae
#168009 syrMik
#156676 tymCupPin

# combine cds seqs per gene
awk '{seq[$1 "\t" $2] = seq[$1 "\t" $2] $3} END {for (i in seq) {print i "\t" seq[i]}}' all_cds.tab > all_cds_final.tab

# convert back to FASTA by second field
mkdir -p unaligned
awk 'length($2) <= 252 {sub("/", "_", $2); print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cds_final.tab

mkdir -p aligned
cd unaligned/

sbatch run_muscle.sh

# if all seqs don't finish with first run, easiest thing to do is to move finished FASTAs to a temp dir and re-run MUSCLE
# ls *.afa > finished.txt
# sed 's/\.afa//g' finished.txt > finished_seq.txt
# mkdir finished_seq/
# cat finished_seq.txt | xargs mv -t finished_seq/
# sbatch run_muscle.sh

mkdir ../aligned
mv *.afa ../aligned
cd ../aligned 

# run HmmCleaner on PRANK alignments
for file in *.afa;
do
  singularity exec --cleanenv /n/singularity_images/informatics/hmmcleaner/hmmcleaner_0.180750.sif HmmCleaner.pl $file
done

