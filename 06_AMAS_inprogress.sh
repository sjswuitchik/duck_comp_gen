module load Anaconda3/5.0.1-fasrc02 
#conda create -n amas python=3.6 amas 
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/alignments/aligned

source activate amas
AMAS.py concat -i batch0*_output/*.aligned.fa -f fasta -d dna -p cnee_partitions.txt -t cnee_concat0.fa -c 12
AMAS.py concat -i batch1*_output/*.aligned.fa -f fasta -d dna -p cnee_partitions1.txt -t cnee_concat1.fa -c 12
AMAS.py concat -i batch2*_output/*.aligned.fa -f fasta -d dna -p cnee_partitions2.txt -t cnee_concat2.fa -c 12
AMAS.py concat -i batch3*_output/*.aligned.fa -f fasta -d dna -p cnee_partitions3.txt -t cnee_concat3.fa -c 12




AMAS.py concat -i cnee_concat*.fa -f fasta -d dna -p cnee_partitions_final.txt -t cnee_concat_final.fa -c 12 

# clean up
cat cnee_concat_final.fa | perl -p -e 's/[?]/-/g' > cnee_concat_gapsFixed.fa
AMAS.py remove -i cnee_concat_gapsFixed.fa -f fasta -d dna -g cnee_final -u fasta -c 12

#fix partition file
cat cnee_partitions*.txt | perl -p -e 's/p\d+[_](zfCNEE\d+)=(\d+)-(\d+)/$1\t$2\t$3/' | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' > cnee_partitions.bed



