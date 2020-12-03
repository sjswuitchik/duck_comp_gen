## Prep, top 3
# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/PhyloAcc

# will run batches of 2000 elements each

# set up batch files
mkdir -p top3_batches
# shuffle lines in random order
shuf -i 0-356467 > top3_batches/full_list 
# make input files
split -d -a 3 -l 2000 top3_batches/full_list top3_batches/batch

# set up parameter and output files
mkdir -p top3_param
mkdir -p top3_outs

for I in $(seq 0 178); 
do
  printf -v BATCH "%03d" $I
  PARTFILE=top3_batches/batch$BATCH
  PREFIX=batch${BATCH}
  cat > top3_param/run$I <<EOF 
PHYTREE_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galloTop3.named.mod
SEG_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/alignments/aligned/galloseq.part.bed
ALIGN_FILE /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/alignments/aligned/galloseq_gapFixed.fa
RESULT_FOLDER top3_outs
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


