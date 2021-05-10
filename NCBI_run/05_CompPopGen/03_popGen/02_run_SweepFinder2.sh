# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/sweepfinder

#conda create -c bioconda -n sweepfinder sweepfinder2 vcftools
source activate sweepfinder
# run `SweepFinder2` to see usage

# copy over BH duck VCF
cp -v ../hetAtr_stiNae_qc/hetAtr.filtered.sorted.vcf.gz* .

# extract allele counts & prep input file for combined allele freq
vcftools --gzvcf hetAtr.filtered.sorted.vcf.gz --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-filtered-all --remove-indels --out hetAtr --counts2
tail -n+2 hetAtr.frq.count | awk -v OFS='\t' '{print $2,$6,$4,"1"}' > hetAtr.combo.sweep
echo -e 'position\tx\tn\tfolded' | cat - hetAtr.combo.sweep > tmp && mv tmp hetAtr.combo.sweep

# create empirical freq spectrum
SweepFinder2 -f hetAtr.combo.sweep hetAtr.spect

## create chromosome-specific allele freq inputs
# delete header row
sed '1d' hetAtr.frq.count > hetAtr.frq 
# split into chr-specific files
awk -F '\t' '{print > $1}' hetAtr.frq
# delete all empty files
find -size 0 -print -delete

## this leaves you with a set of scaffold- and chromosome-specific allele frequency files that can be converted to SweepFinder input

# prep each chromosome allele counts for SweepFinder input as allele freq file
for file in CM*;
do
  tail -n+2 $file | awk -v OFS='\t' '{print $2,$6,$4,"1"}' > $file.sweep
  echo -e 'position\tx\tn\tfolded' | cat - $file.sweep > tmp && mv tmp $file.sweep
done

# run SweepFinder on each main chromosome 



SweepFinder2 -lg 1000 CM021731.1.sweep hetAtr.spect CM021731.1.sweep.out
