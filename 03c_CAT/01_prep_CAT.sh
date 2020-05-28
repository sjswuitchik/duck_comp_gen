##############################
## Gene annotation with CAT ##
##############################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/

# set up CAT
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit.git
pip install -e Comparative-Annotation-Toolkit

# extract only focal and reference species from WGA
cd Comparative-Annotation-Toolkit/
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
hal2maf ../galloanserae.hal galloForCAT.maf --refGenome galGal --noAncestors --noDupes --targetGenomes galGal,hetAtr,netAur,oxyJam,stiNae
maf2hal galloForCAT.maf galloForCAT.hal --refGenome galGal
exit 

export SINGULARITY_IMAGE=/n/singularity_images/informatics/cat/cat:20200116.sif
# run test
luigi --module cat RunCat --hal=test_data/vertebrates.hal --ref-genome=mm10 --workers=10 --config=test_data/test.config --work-dir test_install --out-dir test_install --local-scheduler --augustus --augustus-cgp --augustus-pb --assembly-hub > log.txt







