## FASTQ to VCF with snakemake workflow  
  
This directory contains the files necessary to run an earlier version of the snpArcher workflow (https://github.com/harvardinformatics/snpArcher) of FASTQ to VCF generation. The files and directories contained here include:  
  
* script `01_run_fastq2vcf.sh`, which contains the code for generating a VCF containing black-headed duck reseq data and a male freckled duck aligned to the black-headed duck genome  
* script `02_run_fastq2vcf_reseq.sh`, which contains the code for generating VCFs for each focal species aligned to the respective reference genome
* script `03_vcfQC.sh`, which contains the code for initial QC on the black-headed duck and the freckled duck VCFs, both aligned to the black-headed duck reference genome  
* 












## Notes for running snakemake pipeline ##

### allDucks run ###
10X data needs memory bump for dedup, so for now, aligning all ducks reseq data to hetAtr genome without 10X data from genome asssemblies
- ref: "data/allDucks/genome/hetAtr.fa"
- fastqDir: "data/allDucks/fastqs/" 
- fastq_suffix1: "_1.fastq.gz"
- fastq_suffix2: "_2.fastq.gz"

### species-specific runs ###
All species are aligned to their own genome with the female genome data & the male reseq data (e.g. stiNae female + male aligned to stiNae genome). Black-headed duck (hetAtr) run also includes the male stiNae reseq data for downstream analysis in MK pipeline, with the stiNae male acting as the outgroup.  

NB: 10X data crashes dedup, so for now, not dedup'ing BAMs in hetAtr, netAur, oxyJam, and stiNae run. To do so, in the fastq2bam/01_mappedReads dir:  
- ```cp sample_sorted.bam sample_dedup.bam```  
- ```samtools index sample_dedup.bam```  
- ```mv sample_dedup.bam.bai sample_dedup.bai```  
- ```touch ../02_bamSumstats/sample_dedupMetrics.txt```  

Then re-submit run_fastq2bam_spp.sh to complete the pipeline.
