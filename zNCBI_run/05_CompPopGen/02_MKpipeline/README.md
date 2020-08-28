# Pipeline for comparative genomics 

Automatic pipeline to concatenate and filter VCFs, annotate variants with snpEff, and run SnIPRE, MK tests, and direction of selection calculations. Minor configuration and preprocessing required before running pipeline.sh, outlined below.

Authors: 


Sara Wuitchik (Postdoctoal associate, Boston University & Harvard University; sjswuit@g.harvard.edu)

Allison Shultz (Assistant Curator of Ornithology, LA Natural History Museum; ashultz@nhm.org)

Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)

## Configuration and set up

Project directory should be set up as outlined in directory_tree.pdf so the main pipeline can automatically switch between directories and call files as necessary. 

First, set up a conda environment that will allow access to python and R packages:

```conda create -n mk_v2 python=3.6 anaconda cyvcf2 tqdm bcftools vcftools htslib java-jdk bedtools r-base r-tidyverse r-rjags r-r2jags r-lme4 r-arm```

```source activate mk_v2```

A few variables need to be set before running:

```export R_LIBS_USER=$HOME/path/to/R/packages```

```export PATHW=$HOME/path/to/working/directory```

```export INSHORT=ingroup_spp_name (six letter code)```

```export OUTSHORT=outgroup_spp_name (six letter code)```

```export INLONG=ingroup_spp_name (longform spp name, with leading underscore)```

```export OUTLONG=outgroup_spp_name (longform spp name, with leading underscore)```


Example species names variables: 

```export INSHORT=corCor```

```export OUTSHORT=corMon```

```export INLONG=_Ccornix```

```export OUTLONG=_Cmonedula```

### SnpEff

We use SnpEff (http://snpeff.sourceforge.net/download.html) to build databases and annotate the variants in the VCFs. It should be downloaded in your project directory and set up prior to running the pipeline.

```wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip```

```unzip snpEff_latest_core.zip```

```rm snpEff_latest_core.zip``` 

```rm -r clinEff/```

```cd snpEff/```

```mkdir -p data/$INSHORT/```

Ensure reference sequence (FASTA) and genome annotation (GFF3) are in the appropriate data directory, rename files to sequences.fa and genes.gff, then gzip.

#### Add genome information to config file

Add the following to the snpEff.config file, under the Databases & Genomes - Non-standard Genomes section:

\# Common name genome, Source and Version

$INSHORT.genome : genome_name

For example: 

\# Hooded crow genome, NCBI version 2

corCor.genome : Corvus_cornix_cornix

Once snpEff is ready, export the path to a variable:

```export PATHS=$HOME/path/to/snpEff```

#### Build a snpEff database

From your working directory, run: 

```java -jar $PATHS/snpEff.jar build -gff3 -v $INSHORT```

### In your working directory, you'll need: 

- Missingness data (all_all_missingness_info.txt) for both ingroup and outgroup

- Coverage site data (clean_coverage_sites_merged.bed) for both ingroup and outgroup

- genes.gff (same file that's in the snpEff data directory, gunzipped)

- genenames.py

- gff2bed.awk

- parser_nov.py

- my.jags2.R 

- SnIPRE_source.R

- missingness.R

- mktest.R

- ingroup.vcfs and outgroup.vcfs directories (containing hardfiltered VCFs) 

#### All of the above information can be found in config.sh; once working directory is ready, pipeline.sh can be run from within the directory to output results from SnIPRE, MK tests, and direction of selection calculations
