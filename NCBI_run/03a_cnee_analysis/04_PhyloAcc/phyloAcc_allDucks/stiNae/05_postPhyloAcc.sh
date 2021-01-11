# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/stiNae/

module load bedtools2/2.26.0-fasrc01 Anaconda/5.0.1-fasrc01 R/4.0.2-fasrc01

cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_cnees/galGal6_final_merged_CNEEs_named.bed .

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
# use output from phyloP_cleanup_cnees.R to bin accelerated CNEEs list into windows
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

# sort full list of CNEEs
bedtools sort -i galGal6_final_merged_CNEEs_named.bed > galGal6_final_merged_CNEEs_named_sorted.bed

# annotate slopped BED with accel CNEEs 
bedtools annotate -i galGal.slop.bed -files acc.cnees.final.bed galGal6_final_merged_CNEEs_named_sorted.bed -counts > cnee_gene100kb.bed 
awk '$6 > 0 {print}' cnee_gene100kb.bed | sed 's/ /\t/g' - > cnee_gene100kb.clean.bed

# do permutations
mkdir cnee_perms/
for i in {0001..1000}; 
do
  shuf -n 4376 galGal6_final_merged_CNEEs_named_sorted.bed > cnee_perms/'shuffle'$i.bed
done

# annotate 
bedtools annotate -i cnee_gene100kb.clean.bed -files cnee_perms/*.bed -counts > cnee_perms.counts.bed


### input generation for GO permutations 
# nb: could also write quick script (like replace_chr.pl) to replace the gene symbols with NCBI gene IDs in cnee_perms.counts.bed
# pull out gene symbols & NCBI IDs
mkdir go_perms
cat galGal.gff.bed | python3 gene_ncbi_names.py > go_perms/galGal.ncbigenes.bed
cd go_perms/

# sort 
bedtools sort -i galGal.ncbigenes.bed > galGal.ncbigenes.sorted.bed

# set 100 kb window around genes
bedtools slop -i galGal.ncbigenes.sorted.bed -g ../galGal.chrom.sizes -b 100000 > galGal.ncbislop.bed

# annotate slopped BED with accel CNEEs 
bedtools annotate -i galGal.ncbislop.bed -files ../acc.cnees.final.bed ../galGal6_final_merged_CNEEs_named_sorted.bed -counts > cnee_ncbigene100kb.bed 
awk '$8 > 0 {print}' cnee_ncbigene100kb.bed | sed 's/ //g' - > cnee_ncbigene100kb.clean.bed

# do permutations
mkdir perms/
for i in {0001..1000}; 
do
  shuf -n 4376 ../galGal6_final_merged_CNEEs_named_sorted.bed > perms/'shuffle'$i.bed
done

# annotate 
bedtools annotate -i cnee_ncbigene100kb.clean.bed -files perms/*.bed -counts > cnee_goperms.counts.bed
