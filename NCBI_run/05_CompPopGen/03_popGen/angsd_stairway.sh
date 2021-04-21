# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

source activate vcfqc
# get Stairway
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
unzip stairway_plot_v2.1.1.zip
rm stairway_plot_v2.1.1.zip
mv stairway-plot-v2/ stairway/


# run angsd
sbatch run_angsd.sh


cp ../hetAtr.filtered.vcf .
samtools faidx hetAtr.ncbi.fasta -o hetAtr.ncbi.fasta.fai
vcftools --vcf hetAtr.filtered.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --out hetAtr.stair --recode --recode-INFO-all



