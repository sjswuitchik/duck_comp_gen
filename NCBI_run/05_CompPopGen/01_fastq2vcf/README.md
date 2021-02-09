## Notes for running snakemake pipeline ##

Each run requires changes to the config.yaml that are specific to the run. Below are notes about changes made for each particular run of the pipeline. If a parameter is not listed, it was run with the default setting.

### allDucks run ###
10X data needs memory bump for dedup, so for now, aligning all ducks to hetAtr genome without 10X data from genome asssemblies
- ref: "data/allDucks/genome/hetAtr.fa"
- fastqDir: "data/allDucks/fastqs/" 
- fastq_suffix1: "_1.fastq.gz"
- fastq_suffix2: "_2.fastq.gz"

For bam2vcf
- bamsForGatk: "fastq2bam/01_mappedReads/"
- bam_suffix: "_dedup.bam"

### netAur run ###
NB: 10X data crashes dedup, so for now, not dedup'ing BAMs - manually renaming sorted BAM to dedup'd BAM, indexing, and going on with pipeline
Everything is the same as the above allDucks parameters with the following exception for fastq2bam & intervals
- ref: "data/netAur/genome/netAur.fa"
- fastqDir: "data/netAur/fastqs/" 
- minNmer: 85000
- maxIntervalLen: 210000000 
