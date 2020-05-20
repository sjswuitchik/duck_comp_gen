#################
## PhyloP test ##
#################

module load Anaconda3/2019.10

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/phyloP_test

conda activate phast 
cp -v ../PhyloAcc/input_data/galloseq_gapFixed.fa ../PhyloAcc/input_data/galloseq.part.bed ../PhyloAcc/input_data/galloTop1.named.mod .
# usage: phyloP [OPTIONS] tree.mod [alignment] > out
phyloP --method LRT --mode CON --msa-format FASTA --features galloseq.part.bed galloTop1.named.mod galloseq_gapFixed.fa > cnees_phyloP.out

phyloP --method LRT --mode CON --msa-format FASTA --features galloseq.part.bed galloTop2.named.mod galloseq_gapFixed.fa > cnees_phyloP_top2.out


