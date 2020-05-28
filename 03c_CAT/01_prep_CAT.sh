##############################
## Gene annotation with CAT ##
##############################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/

# set up CAT
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit.git
pip install -e Comparative-Annotation-Toolkit

# extract only focal and reference species from WGA
mkdir -p Comparative-Annotation-Toolkit/input_data
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
hal2maf galloanserae.hal galloForCAT.maf --refGenome galGal --noAncestors --noDupes --targetGenomes galGal,hetAtr,netAur,oxyJam,stiNae
maf2hal galloForCAT.maf galloForCAT.hal --refGenome galGal
exit 

# get input data set up
mv galloForCAT.hal Comparative-Annotation-Toolkit/input_data
cd Comparative-Annotation-Toolkit/input_data/
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/4d_sites/galGal6.filtpy.gff . 



