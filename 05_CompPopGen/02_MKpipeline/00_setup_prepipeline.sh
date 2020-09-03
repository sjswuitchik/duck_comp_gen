# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/MK_pipeline

module load Anaconda3/2019.10
#conda create -c bioconda -c conda-forge -n busco busco=4.0.6 gffread

# only using previously build busco env for the gffread util
conda activate busco

# from /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05b_comppopgen_snakemake/02_MK_pipeline, copy over these files: 
# hetAtr.vcf.gz
# stiNae.vcf.gz
# hetAtr_all_all_missingness_info.txt
# stiNae_all_all_missingness_info.txt

#### still need to figure out what the coverage output from snakemake will look like

cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/output_gtfs/cleaned_reordered_* .
for file in hetAtr netAur oxyJam stiNae;
do
  gffread cleaned_reordered_$file.sorted.gtf -o $file.gff
  
### finish this up after dog walk ...   
