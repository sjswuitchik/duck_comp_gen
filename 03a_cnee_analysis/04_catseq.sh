## Gene concatenation for CNEE analysis ##
## https://github.com/ChrisCreevey/catsequences ## 
## NB: run from git clone  --branch seqname https://github.com/harvardinformatics/catsequences.git until PR is fulfilled ##

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/alignments/aligned

find . -name '*.fa' > list
git clone  --branch seqname https://github.com/harvardinformatics/catsequences.git
cd catsequences/
cc catsequences.c -o catsequences -lm
cd ..
catsequences/catsequences list
mv allseqs.fas galloseq.fa
mv allseqs.partitions.txt galloseq.partitions.txt

# clean up ?s
cat galloseq.fa | perl -p -e 's/[?]/-/g' > galloseq_gapFixed.fa

# need to make part.txt into a bed with CNEE-start-end
sed 's/\.\/batch086_output\///g' galloseq.partitions.txt | awk 'BEGIN{FS="="; OFS="\t"} {split($2,a,"-"); print $1,a[1],a[2]}' | sed 's/;$//' | sed 's/\.aligned\.fa//g' | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' > galloseq.part.bed
