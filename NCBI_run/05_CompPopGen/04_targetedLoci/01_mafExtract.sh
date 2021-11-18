# in /n/holyscratch01/informatics/swuitchik/gallo_align

# convert WGA HAL to MAF
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif
hal2maf gallo_ncbi.hal gallo_ncbi.maf --refGenome galGal --noAncestors --noDupes

# extract one test loci 
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006127.5 --start 65334667 --stop 65337079 > test.maf

# extract all loci of interest 
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006127.5 --start 71189060 --stop 71190152 > ACO1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006127.5 --start 65334667 --stop 65337079 > ALDOB.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006092.5 --start 7824930 --stop 7825369 > ARNTL.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006092.5 --start 3980268 --stop 3980955 > BDNF.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006091.5 --start 65108758 --stop 65109290  > CLOCK.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006106.5 --start 7435280 --stop 7437118 > CLTC.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006102.5 --start 357265 --stop 357902 > CLTC1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006088.5 --start 111075349 --stop 111076534 >CRYAA.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006097.5 --start 9678132 --stop 9678665 > CYP19A1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006115.5 --start 1556621 --stop 1558322 > EER2.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006100.5 --start 18951991 --stop 18952505 > EGR1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006100.5 --start 18950768 --stop 18951969 > EGR1b.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006091.5 --start 20311976 --stop 20312574 > FGB_45.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006091.5 --start 20312599 --stop 20313062 > FGB_56.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006091.5 --start 20313083 --stop 20314577 > FGB_68.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006110.5 --start 270565 --stop 272152 > HMGN2.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006089.5 --start 32590393 --stop 32591849 > HOXA3.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006100.5 --start 17595133 --stop 17596071 > IRF1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006091.5 --start 39543572 --stop 39544188 > IRF2.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006088.5 --start 52006409 --stop 52007352 > MB.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006127.5 --start 66929565 --stop 66930206 > MUSK.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006089.5 --start 139737309 --stop 139738535 > MYC.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006113.5 --start 4027678 --stop 4028426 > NGF.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006088.5 --start 74284706 --stop 74285433 > NTF3.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006093.5 --start 12328860 --stop 12329878 > PCBD1.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006099.5 --start 19917025 --stop 19918579 > RHO.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006127.5 --start 43252453 --stop 43253246 > SPIN.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006090.5 --start 19454853 --stop 19455427 > TGFB2.maf
singularity exec --cleanenv /n/singularity_images/informatics/maftools/maftools:20170913.sif mafExtractor --maf gallo_ncbi.maf --seq galGal.NC_006097.5 --start 4626689 --stop 4627147 > TPM1.maf

# get coordinates and write into BED files
python3 GetCoordinates.py

# copy over reference genomes for FASTA creation
cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/genomes/*.fasta .

# using ref genomes and BEDs, write out FASTAs for each gene 
for file in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
  ./get_fasta.sh $file.ncbi.fasta bed_files/$file.bed
done



