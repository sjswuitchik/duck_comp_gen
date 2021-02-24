## Notes for running snakemake pipeline ##

### allDucks run ###
10X data needs memory bump for dedup, so for now, aligning all ducks reseq data to hetAtr genome without 10X data from genome asssemblies
- ref: "data/allDucks/genome/hetAtr.fa"
- fastqDir: "data/allDucks/fastqs/" 
- fastq_suffix1: "_1.fastq.gz"
- fastq_suffix2: "_2.fastq.gz"


NB: 10X data crashes dedup, so for now, not dedup'ing BAMs in netAur, oxyJam, and stiNae run. To do so, in the fastq2bam/01_mappedReads dir:  
- ```cp sample_sorted.bam to sample_dedup.bam```  
- ```samtools index sample_dedup.bam```  
- ```mv sample_dedup.bam.bai sample_dedup.bai```  
- ```touch ../02_bamSumstats/sample_dedupMetrics.txt```

