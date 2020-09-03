#!/usr/bin/bash

# create an environment that has all the other packages you'll need 
conda create -n mk_v2 python=3.6 anaconda cyvcf2 tqdm bcftools vcftools htslib java-jdk bedtools r-base r-tidyverse r-rjags r-r2jags r-lme4 r-arm
source activate mk_v2

# set up project directory as the structure outlined in directory_tree.pdf 

export R_LIBS_USER=$HOME/path/to/R/packages
export INSHORT=ingroup_spp_name (six letter code)
export OUTSHORT=outgroup_spp_name (six letter code)
export INLONG=ingroup_spp_name (longform spp name, with leading underscore)
export OUTLONG=outgroup_spp_name (longform spp name, with leading underscore)

# for example: 

# export R_LIBS_USER=$HOME/apps/R_3.6.1:$R_LIBS_USER
# export INSHORT=corCor
# export OUTSHORT=corMon
# export INLONG=_Ccornix
# export OUTLONG=_Cmonedula

export PATHW=$HOME/path/to/working/directory

wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip 
rm -r clinEff/
cd snpEff/
mkdir -p data/$INSHORT/
# ensure reference sequence (FASTA) and genome annotation (GFF3) are in the appropriate data directory
# rename reference sequence to sequences.fa and gzip
# rename genome annotation to genes.gff and gzip

# add the following into the snpEff.config file under the Databases & Genomes section: 

# Common name genome, Source and Version
$INSHORT.genome : genome_name

# for example: 

# # Hooded crow genome, NCBI version 2
# corCor.genome : Corvus_cornix_cornix


export PATHS=$HOME/path/to/snpEff

# build database
java -jar $PATHS/snpEff.jar build -gff3 -v $INSHORT

# in working directory, will need: 
# ingroup missingness
# outgroup missingness
# ingroup coverage sites
# outgroup coverage sites
# genes.gff
# genenames.py
# gff2bed.awk
# parser_nov.py
# my.jag2.R
# SnIPRE_source.R
# missingness.R
# mktest.R
# ingroup.vcfs dir (w/ hardfiltered vcfs)
# outgroup.vcfs dir (w/ hardfiltered vcfs)

