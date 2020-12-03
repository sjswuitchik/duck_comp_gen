## Notes for running snakemake pipeline ##

Each run requires changes to the config.yaml that are specific to the run. Below are notes about changes made for each particular run of the pipeline. 

### hetAtr & stiNae run for MK pipline ###
For fastq2bam
- ref: "data/hetAtr/genome/hetAtr.fa"
- fastqDir: "data/hetAtr/fastqs/" 
- fastq_suffix1: "_1.fastq.gz"
- fastq_suffix2: "_2.fastq.gz"

For bam2vcf
- bamsForGatk: "fastq2bam/01_mappedReads/"
- bam_suffix: "_dedup.bam"

### allDucks run ###
10x data needs memory bump for Picard MarkDuplicates, so for now, aligning all ducks to hetAtr genome 
Everything is the same as the above hetAtr parameters with the following exception for the fastq2bam pipeline
- fastqDir: "data/allDucks/fastqs" 

