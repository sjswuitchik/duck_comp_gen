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








# ind01
cat hetAtr_ind01_L1_573_1.fastq.gz hetAtr_ind01_L1_574_1.fastq.gz hetAtr_ind01_L1_575_1.fastq.gz hetAtr_ind01_L1_576_1.fastq.gz hetAtr_ind01_L1_577_1.fastq.gz hetAtr_ind01_L1_578_1.fastq.gz hetAtr_ind01_L1_579_1.fastq.gz hetAtr_ind01_L1_580_1.fastq.gz > hetAtr_ind01_1.fastq.gz

cat hetAtr_ind01_L1_573_2.fastq.gz hetAtr_ind01_L1_574_2.fastq.gz hetAtr_ind01_L1_575_2.fastq.gz hetAtr_ind01_L1_576_2.fastq.gz hetAtr_ind01_L1_577_2.fastq.gz hetAtr_ind01_L1_578_2.fastq.gz hetAtr_ind01_L1_579_2.fastq.gz hetAtr_ind01_L1_580_2.fastq.gz > hetAtr_ind01_2.fastq.gz
