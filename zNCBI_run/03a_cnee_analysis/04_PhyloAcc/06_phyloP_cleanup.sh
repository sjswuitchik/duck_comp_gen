#################
## PhyloP test ##
#################

module load Anaconda3/2019.10

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc

conda activate phast 

# usage: phyloP [OPTIONS] tree.mod [alignment] > out
phyloP --method LRT --mode CON --msa-format FASTA --features input_data/galloseq_ncbi.part.bed input_data/galloTop1.named.mod input_data/galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP_top1.out
phyloP --method LRT --mode CON --msa-format FASTA --features input_data/galloseq_ncbi.part.bed input_data/galloTop2.named.mod input_data/galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP_top2.out
phyloP --method LRT --mode CON --msa-format FASTA --features input_data/galloseq_ncbi.part.bed input_data/galloTop3.named.mod input_data/galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP_top3.out
