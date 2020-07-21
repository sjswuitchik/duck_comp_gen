# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees

module load bedtools2/2.26.0-fasrc01 Anaconda/5.0.1-fasrc01

mkdir postPhyloAcc
cp GCF_000002315.6_GRCg6a_genomic.gff.gz galGal6_final_merged_CNEEs_named.bed postPhyloAcc/
cd postPhyloAcc
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz

# gff to bed
column -s, -t < GCF_000002315.6_GRCg6a_genomic.gff | awk '$3 == "CDS"' > galGal.onlyCDS.gff
awk -f gff2bed.awk galGal.onlyCDS.gff > galGal.gff.bed
tail -n +2 galGal.gff.bed > galGal.clean.gff.bed

# pull out genes
cat galGal.gff.bed | python3 genenames.py > galGal.genes.bed

# sort input
bedtools sort -i galGal.genes.bed > galGal.genes.sorted.bed
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed

# find closest gene to CNEEs
bedtools closest -a galGal6_final_merged_CNEEs_named_sorted.bed -b galGal.genes.sorted.bed | cut -f1,2,3,4,8 | bedtools merge -i - -d -1 -c 4,5 -o distinct> galGal_cnees_genes.bed
