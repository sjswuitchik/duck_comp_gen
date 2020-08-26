# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc

module load bedtools2/2.26.0-fasrc01 Anaconda/5.0.1-fasrc01 R/3.6.3-fasrc01

mkdir postPhyloAcc
cp galGal6_final_merged_CNEEs_named.bed postPhyloAcc/
cd postPhyloAcc/

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



