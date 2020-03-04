#### HAL Tools summaries and MAF conversion for CESAR ####
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus on bioinf01

git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit.git
cp ../ducks_cactus/galloanserae.hal Comparative-Annotation-Toolkit/
cd Comparative-Annotation-Toolkit
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
halValidate duck_data/galloanserae.hal
hal2maf data/galloanserae.hal data/galloanserae.maf
halStats data/galloanserae.hal > data/halStats.out
halSummarizeMutations data/galloanserae.hal > data/halSumMut.out
exit

#### Building a set of consensus CNEEs from literature #### 
# get data from /n/holylfs/LABS/informatics/tsackton/broodParasites/DATA/02_CNEEs/concat_cnees to build consensus CNEE set
# notes on data downloads can be found at https://github.com/tsackton/brood-parasite-genomics/tree/master/03_CNEEs/assemble_ces.txt
cp FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol.bed galGal4_phastCons_SP_conserved_elements_LoweEtAlMBE.bed galGal4_phastCons_TP_conserved_elements_SacktonEtAl.bed hg38_phastCons_SP_conserved_elements_UCSC.bed replace_chrs.pl /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees
cd /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees
module load perl/5.26.1-fasrc01 python/2.7.14-fasrc02 bedtools2/2.26.0-fasrc01 Anaconda3/5.0.1-fasrc02 

# make chr keys 
# ficAlb
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/815/GCF_000247815.1_FicAlb1.5/GCF_000247815.1_FicAlb1.5_assembly_report.txt
sed 's/\r$//g' GCF_000247815.1_FicAlb1.5_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > ficAlb_chr_key
# galGal4
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.3_Gallus_gallus-4.0/GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt
sed 's/\r$//g' GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > galGal4_chr_key
# galGal5
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.4_Gallus_gallus-5.0/GCF_000002315.4_Gallus_gallus-5.0_assembly_report.txt
sed 's/\r$//g' GCF_000002315.4_Gallus_gallus-5.0_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > galGal5_chr_key
# galGal6
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt
sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > galGal6_chr_key

# register for liftOver product on UCSC website (free for academic use)
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
# get liftOver chains
# galGal 4 -> galGal6
wget http://hgdownload.soe.ucsc.edu/goldenPath/galGal4/liftOver/galGal4ToGalGal6.over.chain.gz
gunzip galGal4ToGalGal6.over.chain.gz
# hg38 -> galGal6 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToGalGal6.over.chain.gz
gunzip hg38ToGalGal6.over.chain.gz
# galGal5 -> galGal6
wget http://hgdownload.soe.ucsc.edu/goldenPath/galGal5/liftOver/galGal5ToGalGal6.over.chain.gz
gunzip galGal5ToGalGal6.over.chain.gz
# liftOver <CNEEs .bed file> <liftover_chain> <out.bed> <unmapped.out.bed>

# liftOver Lowe et al galGal4 -> galGal6
./liftOver galGal4_phastCons_SP_conserved_elements_LoweEtAlMBE.bed galGal4ToGalGal6.over.chain galGal6_phastcons_SP_conserved_elements_LoweEtAlMBE.bed LoweEtAlUnmapped.bed
awk '{print $4, $3}' galGal6_chr_key > acckey
./replace_chrs.pl acckey galGal6_phastcons_SP_conserved_elements_LoweEtAlMBE.bed > galGal6_phastcons_SP_conserved_elements_LoweEtAlMBE_NCBI.bed

# liftOver Sackton et al galGal4 -> galGal6
cut -f3,4 galGal4_chr_kay > acckey
./replace_chrs.pl acckey galGal4_phastCons_TP_conserved_elements_SacktonEtAl.bed > galGal4_phastCons_TP_conserved_elements_SacktonEtAl_UCSC.bed 
./liftOver galGal4_phastCons_TP_conserved_elements_SacktonEtAl_UCSC.bed galGal4ToGalGal6.over.chain galGal6_phastcons_TP_conserved_elements_SacktonEtAl.bed SacktonEtAlUnmapped.bed
awk '{print $4, $3}' galGal6_chr_key > acckey
./replace_chrs.pl acckey galGal6_phastcons_TP_conserved_elements_SacktonEtAl.bed > galGal6_phastcons_TP_conserved_elements_SacktonEtAl_NCBI.bed

# liftOver UCSC 100-way vertebrate CNEEs, referenced on Homo sapiens, to galGal6
./liftOver hg38_phastCons_SP_conserved_elements_UCSC.bed hg38ToGalGal6.over.chain galGal6_phastcons_SP_conserved_elements_UCSC.bed UCSC_hsap_unmapped.bed
awk '{print $4, $3}' galGal6_chr_key > acckey
./replace_chrs.pl acckey galGal6_phastcons_SP_conserved_elements_UCSC.bed > galGal6_phastcons_SP_conserved_elements_UCSC_NCBI.bed

# replace chr in Craig et al set (eg/ chr 1 -> NW_004775899.1) and use Tim's alignment to liftover to galGal5, then do liftOver to galGal6

cut -f2,3 ficAlb_chr_key > acckey
./replace_chrs.pl acckey FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol.bed > FicAlb.part1.bed
cut -f1,3 ficAlb_chr_key > acckey
./replace_chrs.pl acckey FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol.bed > FicAlb.part2.bed
cat FicAlb.part1.bed FicAlb.part2.bed > FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol_NCBI.bed
cp FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol_NCBI.bed ../Comparative-Annotation-Toolkit
cd ../Comparative-Annotation-Toolkit
cut -f1,2,3,4 FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol_NCBI.bed > FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol_NCBI_cut.bed
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
halLiftover --noDupes data/broodParaAlign.hal ficAlb FicAlb1.5_phastCons_TP_conserved_elements_CraigEtAlMolEcol_NCBI_cut.bed galGal galGal5_Craig.bed 2> craig_liftover.log
exit
mv galGal5_Craig.bed ../cnees
awk '{print $3, $1}' galGal5_chr_key > acckey
./replace_chrs.pl acckey galGal5_Craig.bed > galGal5_Craig_chr.bed
./liftOver galGal5_Craig_chr.bed galGal5ToGalGal6.over.chain galGal6_Craig.bed Craig_unmapped.bed

# merge within 5 bp
bedtools sort -i galGal6_phastcons_SP_conserved_elements_LoweEtAlMBE_NCBI.bed | bedtools merge -i - -d 5 > galGal6_Lowe_merged.bed
bedtools sort -i galGal6_phastcons_TP_conserved_elements_SacktonEtAl_NCBI.bed | bedtools merge -i - -d 5 > galGal6_Sackton_merged.bed
bedtools sort -i galGal6_phastcons_SP_conserved_elements_UCSC_NCBI.bed | bedtools merge -i - -d 5 > galGal6_UCSC_merged.bed
bedtools sort -i galGal6_Craig.bed | bedtools merge -i - -d 5 > galGal6_Craig_merged.bed 

# get ce lengths 
awk '{print $3-$2, "\tLowe"}' galGal6_Lowe_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tSackton"}' galGal6_Sackton_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tUCSC"}' galGal6_UCSC_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tCraig"}' galGal6_Craig_merged.bed  >> ce.lengths

cat galGal6_*_merged.bed | bedtools sort -i - | bedtools merge -i - | awk '{if ($3-$2 >= 50) print $0}' > galGal6_all_merged.bed

awk '{print $3-$2}' galGal6_all_merged.bed > galGal6_allce.lengths

# get exons
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz
awk 'BEGIN{OFS="\t";} {if ($3 ~ /exon/) print $1, $4-1, $5}' GCF_000002315.6_GRCg6a_genomic.gff | bedtools sort -i - | bedtools merge -i - > galGal6.exon.bed

# get CNEEs
bedtools intersect -v -a galGal6_all_merged.bed -b galGal6.exon.bed > galGal6_final_merged_CNEEs.bed

# add names
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="zfCNEE"NR; print}' galGal6_final_merged_CNEEs.bed > galGal6_final_merged_CNEEs_named.bed

