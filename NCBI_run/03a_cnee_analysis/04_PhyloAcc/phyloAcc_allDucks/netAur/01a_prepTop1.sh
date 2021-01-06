# Prep top 1 for PhyloAcc
# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/netAur

module load gcc/7.1.0-fasrc01 armadillo/7.800.2-fasrc02 gsl/2.4-fasrc01
git clone https://github.com/xyz111131/PhyloAcc.git
cd PhyloAcc/
make

# bring over data for all topologies
mkdir input_data/
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/stiNae/PhyloAcc/input_data/* input_data/

# run batches of 2000 elements each
mkdir -p top1_batches
# shuffle lines in random order with 0-(n-1) from elements
shuf -i 0-375589 > top1_batches/full_list
# make input files
split -d -a 3 -l 2000 top1_batches/full_list top1_batches/batch

# set up parameter and output directories
mkdir -p top1_param
mkdir -p top1_outs

for I in $(seq 0 187); # number of batches generated by line 20
do
  printf -v BATCH "%03d" $I
  PARTFILE=top1_batches/batch$BATCH
  PREFIX=batch${BATCH}
  cat > top1_param/run$I <<EOF 
PHYTREE_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/netAur/PhyloAcc/input_data/galloTop1.named.mod
SEG_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/netAur/PhyloAcc/input_data/galloseq_ncbi.part.bed
ALIGN_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/netAur/PhyloAcc/input_data/galloseq_ncbi_gapFixed.fa
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
TARGETSPECIES netAur
CONSERVE stiNae;hetAtr;oxyJam;anaPla;braCan;ansInd;ansCyg;ansBra;numMel;colVir;tymCupPin;syrMik;galGal;cotJap
GAPCHAR -
NUM_THREAD 8
VERBOSE 0
CONSTOMIS 0.01
GAP_PROP 0.8
TRIM_GAP_PERCENT 0.8
EOF
done
