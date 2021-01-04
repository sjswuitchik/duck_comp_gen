########################
## BUSCO for CESAR QC ##
########################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/

module load Anaconda3/2019.10
conda activate busco

# copy GFFs from compAug output dir
mkdir -p busco/input_data

for file in hetAtr netAur oxyJam stiNae;
do
  cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/augCGP_rnahints/joined_pred/$file.gff busco/
  gffread busco/$file.gff -g genomes/$file.ncbi.fasta -y bucso/input_data/$file.translated.fa -S 
done

