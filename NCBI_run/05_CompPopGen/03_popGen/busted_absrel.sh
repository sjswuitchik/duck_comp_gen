# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy bedtools bioawk
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

# for NCBI annotations, create CDS-only GFFs
for file in galGal ansCyg cotJap numMel anaPla;
do
  column -s, -t < $file.gff | awk '$3 == "CDS"' > $file.cds.gff
  echo -e '##gff-version 3' | cat - $file.cds.gff > tmp && mv tmp $file.cds.gff
done

# extract nucleotide sequences associated with CDS
cd ..
for file in ansBra ansInd braCan colVir hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  bedtools getfasta -fi fastas/$file.ncbi.fasta -bed gffs/$file.gff -name -s > $file.out.fa
done

for file in galGal ansCyg cotJap numMel anaPla;
do
  bedtools getfasta -fi fastas/$file.ncbi.fasta -bed gffs/$file.cds.gff -name -s > $file.out.fa
done

# fix up with bioawk
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
  bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$file"'", $name, $seq}' $file.out.fa >> all_cds.fa
done

# check 
cut -f1,1 all_cds.fa | sort | uniq -c 
#607700 anaPla
#225126 ansBra
#391902 ansCyg
#223647 ansInd
#229390 braCan
#187930 colVir
#566194 cotJap
#669824 galGal
#231898 hetAtr
#231254 netAur
#569903 numMel
#240390 oxyJam
#216192 stiNae
#225650 syrMik
#205651 tymCupPin

# convert back to a FASTA 
mkdir unaligned
awk '{print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cds.fa 
