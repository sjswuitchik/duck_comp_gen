#####################################
## Concatenation for CNEE analysis ## 
#####################################

## NB: run from git clone  --branch seqname https://github.com/harvardinformatics/catsequences.git until PR is fulfilled

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/alignments/aligned

# every dir/file in aligned needs to be target ie/ no non-target files in dir for list
mkdir slurm_outs
mv slurm-* run_mafft.sh slurm_outs/
mv slurm_outs/ ..

# set up and run catseq
find . -name '*.fa' > list
git clone  --branch seqname https://github.com/harvardinformatics/catsequences.git
cd catsequences/
cc catsequences.c -o catsequences -lm
cd ..
catsequences/catsequences list
mv allseqs.fas galloseq_ncbi.fa
mv allseqs.partitions.txt galloseq_ncbi.partitions.txt

# clean up ?s
cat galloseq.fa | perl -p -e 's/[?]/-/g' > galloseq_ncbi_gapFixed.fa

# need to make part.txt into a bed with CNEE-start-end
sed 's/\.\/batch.*_output\///g' galloseq_ncbi.partitions.txt | awk 'BEGIN{FS="="; OFS="\t"} {split($2,a,"-"); print $1,a[1],a[2]}' | sed 's/;$//' | sed 's/\.aligned\.fa//g' | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' > galloseq_ncbi.part.bed

mv ../slurm_outs/ .
