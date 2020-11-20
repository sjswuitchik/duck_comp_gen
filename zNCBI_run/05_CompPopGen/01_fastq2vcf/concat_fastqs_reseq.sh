#!/usr/bin/bash

# netAur
cd netAur/fastqs/
cat netAur_female_L001_R1_001.fastq.gz netAur_female_L004_R1_001.fastq.gz > netAur_female_R1.fastq.gz
cat netAur_female_L001_R2_001.fastq.gz netAur_female_L004_R2_001.fastq.gz > netAur_female_R2.fastq.gz

# oxyJam
cd ../../oxyJam/fastqs
cat oxyJam_female_L001_R1_001.fastq.gz oxyJam_female_L004_R1_001.fastq.gz > oxyJam_female_R1.fastq.gz
cat oxyJam_female_L001_R2_001.fastq.gz oxyJam_female_L004_R2_001.fastq.gz > oxyJam_female_R2.fastq.gz

# stiNae
cd ../../stiNae/fastqs
cat _______ > stiNae_female_R1.fastq.gz
cat _______ > stiNae_female_R2.fastq.gz

# hetAtr
cd ../../hetAtr_v2/fastqs
cat DGAB-CNB0004-CN4-lib1_S1_L004_R1_001.fastq.gz DGAB-CNB0004-CN4-lib1_S3_L001_R1_001.fastq.gz > hetAtr_female_R1.fastq.gz
cat DGAB-CNB0004-CN4-lib1_S1_L004_R2_001.fastq.gz DGAB-CNB0004-CN4-lib1_S3_L001_R2_001.fastq.gz > hetAtr_female_R2.fastq.gz
