## Prep, top 1
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis

module load gcc/7.1.0-fasrc01 armadillo/7.800.2-fasrc02 gsl/2.4-fasrc01
git clone https://github.com/xyz111131/PhyloAcc.git
cd PhyloAcc/
make

mkdir -p PhyloAcc/input_data
cd PhyloAcc/input_data
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/4d_neut_models/galloTop*.named.mod /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/alignments/aligned/galloseq_ncbi.part.bed /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/alignments/aligned/galloseq_ncbi_gapFixed.fa .
cd ..

# run batches of 2000 elements each
mkdir -p top1_batches
# get number of elements
wc -l input_data/galloseq_ncbi.part.bed
# shuffle lines in random order with 0-(n-1) from wc -l above 
shuf -i 0-375589 > top1_batches/full_list
# make input files
split -d -a 3 -l 2000 top1_batches/full_list top1_batches/batch

# set up parameter and output files
mkdir -p top1_param
mkdir -p top1_outs

for I in $(seq 0 187); # number of batches generated by line 20
do
  printf -v BATCH "%03d" $I
  PARTFILE=top1_batches/batch$BATCH
  PREFIX=batch${BATCH}
  cat > top1_param/run$I <<EOF 
PHYTREE_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc/input_data/galloTop1.named.mod
SEG_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc/input_data/galloseq_ncbi.part.bed
ALIGN_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc/input_data/galloseq_ncbi_gapFixed.fa
RESULT_FOLDER top1_outs
PREFIX $PREFIX
ID_FILE $PARTFILE
CHAIN 1
BURNIN 1000
MCMC 4000
CONSERVE_PRIOR_A 5
CONSERVE_PRIOR_B 0.04
ACCE_PRIOR_A 10
ACCE_PRIOR_B 0.2
HYPER_GRATE_A 3
HYPER_GRATE_B 1
TARGETSPECIES hetAtr
CONSERVE netAur;stiNae;oxyJam;anaPla;braCan;ansInd;ansCyg;ansBra;numMel;colVir;tymCupPin;syrMik;galGal;cotJap
GAPCHAR -
NUM_THREAD 8
VERBOSE 0
CONSTOMIS 0.01
GAP_PROP 0.8
TRIM_GAP_PERCENT 0.8
EOF
done
