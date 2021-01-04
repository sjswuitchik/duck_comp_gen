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

# get config.ini from BUSCO git and set up 
wget https://gitlab.com/ezlab/busco/-/raw/master/config/config.ini
# copy each path for respective command, being sure to not include the command in the path - uncomment config lines as necessary (config.ini is in NCBI_run/03d_CompAug)
whereis tblastn
whereis makeblastdb
whereis augustus #this path can be pasted for remaining commands & Augustus scripts in config - nb: need to add scripts/ after path for .pl scripts

# run BUSCO
export BUSCO_CONFIG_FILE=/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03d_CompAug/busco/config.ini
mkdir busco_outs

for file in hetAtr netAur stiNae oxyJam;
do
  busco -m protein -i input_data/$file.translated.fa -o $file -l aves_odb10
done

