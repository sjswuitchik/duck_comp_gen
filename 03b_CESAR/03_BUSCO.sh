########################
## BUSCO for CESAR QC ##
########################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos/busco

#conda create -c bioconda -c conda-forge -n busco busco=4.0.6
conda activate busco

# bring over input data 
mkdir input_data
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos/from_cesar/input_data/* input_data/

# set up config
wget https://gitlab.com/ezlab/busco/-/raw/master/config/config.ini?inline=false
# copy each of these paths into the proper paths for each command, and uncomment lines as necessary (see 03b_CESAR/config.ini for parameters used)
whereis tblastn
whereis makeblastdb
whereis augustus #this path can be pasted for remaining scripts in config aside from sepp

# run BUSCO
export BUSCO_CONFIG_FILE=/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/orthos/busco/config.ini
mkdir busco_outs
for file in hetAtr netAur stiNae oxyJam;
do
  busco -m protein -i input_data/$file.translated.fa -o $file -l aves_odb10
done
