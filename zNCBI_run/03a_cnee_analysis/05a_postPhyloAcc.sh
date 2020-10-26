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



### input generation for the assessment of genes with evidence for excess of nearby accelerated CNEEs

# get galGal6 annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# gff to bed
column -s, -t < GCF_000002315.6_GRCg6a_genomic.gff | awk '$3 == "gene"' > galGal.genes.gff
awk -f gff2bed.awk galGal.genes.gff > galGal.gff.bed

# pull out genes
cat galGal.gff.bed | python3 genenames.py > galGal.genes.bed

# sort 
bedtools sort -i galGal.genes.bed > galGal.genes.sorted.bed

# set 100 kb window around genes
bedtools slop -i galGal.genes.sorted.bed -g galGal.chrom.sizes -b 100000 > galGal.slop.bed

# annotate slopped BED with accel CNEEs 
bedtools annotate -i galGal.slop.bed -files acc.cnees.final.bed galGal6_final_merged_CNEEs_named_sorted.bed -counts > cnee_gene100kb.bed 
awk '$6 > 0 {print}' cnee_gene100kb.bed | sed 's/ /\t/g' - > cnee_gene100kb.clean.bed

# do permutations
mkdir cnee_perms/
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed
for i in {0001..1000}; 
do
  shuf -n 294 galGal6_final_merged_CNEEs_named_sorted.bed > cnee_perms/'shuffle'$i.bed
done

# annotate 
bedtools annotate -i cnee_gene100kb.clean.bed -files cnee_perms/*.bed -counts > cnee_perms.counts.bed



# input generation for cluster profiler (GO enrichment permutations) 
awk '$6 > 0 {print}' cnee_gene100kb.bed | awk '$5 != 0 {$5 = 1} {print}' - | sed 's/ /\t/g' - > cnee_gene100kb.GO.bed

mkdir go_perms/
for i in {0001..1000}; 
do
  shuf -n 294 galGal6_final_merged_CNEEs_named_sorted.bed > go_perms/'go_shuffle'$i.bed
done






# find closest gene to CNEEs
bedtools closest -a galGal6_final_merged_CNEEs_named_sorted.bed -b galGal.genes.sorted.bed | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > galGal_cnees_genes.bed

