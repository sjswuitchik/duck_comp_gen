# in /n/holyscratch01/informatics/swuitchik/ducks/compGen

source activate vcfqc
# get Stairway
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
unzip stairway_plot_v2.1.1.zip
rm stairway_plot_v2.1.1.zip
mv stairway-plot-v2/ stairway/

# run angsd
cp /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/03_CompAugAnnotation/genomes/hetAtr.ncbi.fasta .
samtools faidx hetAtr.ncbi.fasta -o hetAtr.ncbi.fasta.fai

sbatch run_angsd.sh

# run Stairbuilder
mkdir plots
cd stairway_plot_v2.1.1
java -cp stairway_plot_es Stairbuilder hetAtr.blueprint
chmod +x hetAtr.blueprint.sh
./hetAtr.blueprint.sh



vcftools --vcf hetAtr.filtered.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --remove-indels --out hetAtr.stair --recode --recode-INFO-all



