## in 

git clone https://github.com/veg/hyphy-analyses.git

conda activate align

for file in og_fastas/*.fa;
do
  hyphy hyphy-analyses/codon-msa/pre-msa.bf --input $file
  muscle -in ${file}_protein.fas -out ${file}_protein.msa
  hyphy hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${file}_protein.msa --nucleotide-sequences ${file}_nuc.fas --output ${file}_codon.msa --compress No
done
