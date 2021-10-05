# in /n/holyscratch01/informatics/swuitchik/gallo_align

# convert WGA HAL to MAF
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif
hal2maf gallo_ncbi.hal gallo_ncbi.maf --refGenome galGal --noAncestors --noDupes

singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq NC_006127.5 --start 65334667 --stop 65337079
