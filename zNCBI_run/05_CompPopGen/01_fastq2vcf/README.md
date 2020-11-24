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

### hetAtr_v2 run ###
Everything is the same as the above hetAtr parameters with the following exception for the fastq2bam pipeline
- fastqDir: "data/hetAtr_v2/fastqs" 
- 10x data seems to require a memory bump for Picard MarkDuplicates and GATK, so running with mem increase to 40000 for dedup (currently testing mem limits, updated Nov 24)
