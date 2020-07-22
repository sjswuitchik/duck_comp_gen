# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc

module load bedtools2/2.26.0-fasrc01 Anaconda/5.0.1-fasrc01 R/3.6.3-fasrc01

mkdir postPhyloAcc
cp galGal6_final_merged_CNEEs_named.bed postPhyloAcc/
cd postPhyloAcc/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# gff to bed
column -s, -t < GCF_000002315.6_GRCg6a_genomic.gff | awk '$3 == "CDS"' > galGal.onlyCDS.gff
awk -f gff2bed.awk galGal.onlyCDS.gff > galGal.gff.bed

# pull out genes
cat galGal.gff.bed | python3 genenames.py > galGal.genes.bed

# sort input
bedtools sort -i galGal.genes.bed > galGal.genes.sorted.bed
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed

# find closest gene to CNEEs
bedtools closest -a galGal6_final_merged_CNEEs_named_sorted.bed -b galGal.genes.sorted.bed | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct > galGal_cnees_genes.bed

# copy accelerated CNEEs BED from R output to intersect with galGal_cnees_genes.bed
bedtools intersect -a acc_cnees.bed -b galGal_cnees_genes.bed -wb | cut -f1,2,3,4,9 > acc_cnees_genes.bed





