## PhyloAcc for CNEE analysis ## 
## https://github.com/xyz111131/PhyloAcc ## 

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cness
module load gcc/7.1.0-fasrc01 armadillo/7.800.2-fasrc02 gsl/2.4-fasrc01
git clone https://github.com/xyz111131/PhyloAcc.git
cd PhyloAcc
nano ~/.bashrc
### add these to your profile
export LD_LIBRARY_PATH=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/phyloAcc:$LD_LIBRARY_PATH
export LD_RUN_PATH=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cness/phyloAcc:$LD_RUN_PATH
###
source ~/.bashrc
### edit the following in the Makefile
whereis gsl
GSL_INCLUDE=/usr/include/gsl 
GSL_LIB=/usr/share/man/man3/gsl.3.gz
###
make 
sudo make install
# enter password for sudo 

mkdir gallo_results_top1
mkdir gallo_results_top2
mkdir gallo_results_top3

sbatch top1_batch.sh 
sbatch top2_batch.sh 
sbatch top3_batch.sh 

