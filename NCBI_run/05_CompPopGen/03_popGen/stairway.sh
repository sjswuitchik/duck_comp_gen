# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

module load jdk/10.0.1-fasrc01
# get Stairway
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
unzip stairway_plot_v2.1.1.zip
rm stairway_plot_v2.1.1.zip
mv stairway-plot-v2/ stairway/

# filter VCF for allele counts
cd stairway
cp ../hetAtr_stiNae_qc/hetAtr.filtered.vcf .
vcftools --vcf hetAtr.filtered.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --out hetAtr --counts2
# edit frq.counts header to be CHROM  POS N_ALLELES N_CHR MAJOR MINOR
Rscript calc_sfs.R

# run Stairbuilder
cd stairway_plot_v2.1.1
java -cp stairway_plot_es Stairbuilder hetAtr.blueprint
chmod +x hetAtr.blueprint.sh
./hetAtr.blueprint.sh





