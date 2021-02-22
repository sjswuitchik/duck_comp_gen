## Notes for running snakemake pipeline ##

Each run requires changes to the config.yaml that are specific to the run. Below are notes about changes made for each particular run of the pipeline. If a parameter is not listed, it was run with the default setting.

### allDucks run ###
10X data needs memory bump for dedup, so for now, aligning all ducks reseq data to hetAtr genome without 10X data from genome asssemblies
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

### hetAtr + stiNae run for MK pipeline ###  
data: non-dedup'd 10X BAMs, hetAtr reseq data, hetAtr 10X genome data, stiNae male reseq data  
- ref: "data/hetSti/genome/hetAtr.fa"  
- fastqDir: "data/hetSti/fastqs/"

### stiNae run ###
data: non-dedup'd BAMs, female 10X genome data, male reseq data  
- ref: "data/stiNae/genome/stiNae.fa"
- fastqDir: "data/stiNae/fastqs/"  

### oxyJam run ###
data: non-dedup'd BAMs, female 10X genome data, male reseq data  
- ref: "data/oxyJam/genome/oxyJam.fa"
- fastqDir: "data/oxyJam/fastqs/"
