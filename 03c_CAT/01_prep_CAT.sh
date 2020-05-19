##############################
## Gene annotation with CAT ##
##############################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/

singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif

# extract only focal and reference species from WGA
hal2maf galloanserae.hal galloForCAT.maf --refGenome galGal --noAncestors --noDupes --targetGenomes galGal,hetAtr,netAur,oxyJam,stiNae

