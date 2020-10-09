# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc

module load bedtools2/2.26.0-fasrc01 Anaconda/5.0.1-fasrc01 R/4.0.2-fasrc01

mkdir postPhyloAcc
cp galGal6_final_merged_CNEEs_named.bed postPhyloAcc/
cd postPhyloAcc/

### spatial enrichment analyses input generation

# sort full final list
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed

# make galGal6 chrom.sizes file for bedtools makewindows with faToTwoBit and twoBitInfo
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz
gunzip GCF_000002315.6_GRCg6a_genomic.fna.gz

./faToTwoBit GCF_000002315.6_GRCg6a_genomic.fna galGal.fa.2bit
./twoBitInfo galGal.fa.2bit stdout | sort -k2rn > galGal.chrom.sizes

# create 100kb windows with 50kb slide
bedtools makewindows -g galGal.chrom.sizes -w 100000 -s 50000 > galGal.windows.bed
# bin total CNEEs list into windows
bedtools intersect -a galGal.windows.bed -b galGal6_final_merged_CNEEs_named_sorted.bed -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > window.cnees.bed
# use output from 04_PhyloAcc/07_phyloP_cleanup_cnees.R to bin accelerated CNEEs list into windows
bedtools intersect -a galGal.windows.bed -b acc.cnees.final.bed -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > window.acc.cnees.bed

# repeat above steps in 5Mb windows, no slide
bedtools makewindows -g galGal.chrom.sizes -w 5000000 > galGal.5Mbwindows.bed
bedtools intersect -a galGal.5Mbwindows.bed -b galGal6_final_merged_CNEEs_named_sorted.bed -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > 5Mbwindow.cnees.bed
bedtools intersect -a galGal.5Mbwindows.bed -b acc.cnees.final.bed -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > 5Mbwindow.acc.cnees.bed


### cluster profiler input generation

# get galGal6 genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# gff to bed
column -s, -t < GCF_000002315.6_GRCg6a_genomic.gff | awk '$3 == "CDS"' > galGal.onlyCDS.gff
awk -f gff2bed.awk galGal.onlyCDS.gff > galGal.gff.bed

# pull out genes
cat galGal.gff.bed | python3 genenames.py > galGal.genes.bed
awk '{ if (NF == 4) { print } }' galGal.genes.bed > galGal.genes.clean.bed

export R_LIBS_USER=$HOME/apps/R_4.0.2
Rscript --slave --vanilla tidygenes.R 'galGal.NCBIgenes.clean.bed' >std.Rout 2> std.Rerr

# sort 
bedtools sort -i galGal.tidygenes.bed > galGal.genes.sorted.bed

# set 100 kb window around genes
bedtools slop -i galGal.genes.sorted.bed -g galGal.chrom.sizes -b 100000 > galGal.slop.bed

# annotate slopped BED with accel CNEEs 
bedtools annotate -i galGal.slop.bed -files acc_cnees.bed galGal6_final_merged_CNEEs_named_sorted.bed -counts > cnee_gene100kb.bed 
awk '{FS = OFS = '\t'} $6 > 0 {print}' cnee_gene100kb.bed > cnee_gene100kb.clean.bed

# replace '-' in gene names (because FS were being weird) and non-zero accel counts with '1's 
sed -i 's/' ## need to work on this 
awk '{FS = OFS = '\t'} $5 != 0 {$5 = 1} {print}' cnee_gene100kb.clean.bed > cnee_gene100kb.rep.bed

# do permutations
mkdir cnee_perms/
for i in {0001..1000}; 
do
  shuf -n 294 cnee_gene100kb.rep.bed > cnee_perms/'shuffle'$i.bed
done

cd cnee_perms
for file in *.bed;
do
  awk '{FS = OFS = '\t'} (NF == 6) {print}' $file > $file.clean
done
rm *.bed 
rename -v 'bed.clean' 'bed' cnee_perms/*.clean

# annotate perms together
mv cnee_perms/shuffle0001.bed .
cd cnee_perms
ls *.bed > shuffle
mv shuffle ..
cd ..
bedtools annotate -i shuffle0001.bed -files shuffle > cnee_perms.bed

# find closest gene to CNEEs
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed
bedtools closest -a galGal6_final_merged_CNEEs_named_sorted.bed -b galGal.genes.sorted.bed | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > galGal_cnees_genes.bed

