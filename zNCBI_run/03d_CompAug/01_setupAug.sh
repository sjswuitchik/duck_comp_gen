## in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/

module load Anaconda/5.0.1-fasrc01
#conda create -n compAug dockermake boost gsl lpsolve55 pymysql bamtools bcftools htslib samtools
source activate compAug

# install Comp Aug
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus/
make


