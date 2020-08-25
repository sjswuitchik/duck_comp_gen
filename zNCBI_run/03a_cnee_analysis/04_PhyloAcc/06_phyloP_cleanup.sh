#################
## PhyloP test ##
#################

module load Anaconda3/2019.10

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc

conda activate phast 
cp -v input_data/galloseq_gapFixed.fa ../PhyloAcc/input_data/galloseq.part.bed ../PhyloAcc/input_data/galloTop*.named.mod .

# usage: phyloP [OPTIONS] tree.mod [alignment] > out
phyloP --method LRT --mode CON --msa-format FASTA --features galloseq_ncbi.part.bed galloTop1.named.mod galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP.out
phyloP --method LRT --mode CON --msa-format FASTA --features galloseq_ncbi.part.bed galloTop2.named.mod galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP_top2.out
phyloP --method LRT --mode CON --msa-format FASTA --features galloseq_ncbi.part.bed galloTop3.named.mod galloseq_ncbi_gapFixed.fa > cnees_ncbi_phyloP_top3.out
