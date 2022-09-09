## FASTQ to VCF with snakemake workflow  
  
This directory contains the files necessary to run an earlier version of the snpArcher workflow (https://github.com/harvardinformatics/snpArcher) of FASTQ to VCF generation. The files and directories contained here include:  
  
* script `01_run_fastq2vcf.sh`, which contains the code for generating a VCF containing black-headed duck reseq data and a male freckled duck aligned to the black-headed duck genome  
* script `02_run_fastq2vcf_reseq.sh`, which contains the code for generating VCFs for each focal species aligned to the respective reference genome
* script `03_vcfQC.sh`, which contains the code for initial QC on the black-headed duck and the freckled duck VCFs, both aligned to the black-headed duck reference genome  
* `outputs` directory, which contains the outputs from `03_vcfQC.sh`
* `required_files` directory, which contains all snakefiles, run files, and configs for various FASTQ to VCF runs
* `subscripts` directory, which contains any scripts used in the main scripts detailed above
  
  
  
  
#### Species-specific runs ####
All species are aligned to their own genome with the female genome data & the male reseq data (e.g. stiNae female + male aligned to stiNae genome). Black-headed duck (hetAtr) run also includes the male stiNae reseq data for downstream analysis in MK pipeline, with the stiNae male acting as the outgroup.
