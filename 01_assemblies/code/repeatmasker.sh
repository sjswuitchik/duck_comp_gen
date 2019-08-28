module load RepeatMasker/4.0.5-fasrc05

nohup RepeatMasker -species chicken -gff -dir hetAtr_RM/ hetAtr.1.fasta &> nohup_hetAtr_RM.out& 

nohup RepeatMasker -species chicken -gff -dir netAur_RM/ netAur.1.fasta &> nohup_netAur_RM.out& 

nohup RepeatMasker -species chicken -gff -dir oxyJam_RM/ oxyJam.1.fasta.gz &> nohup_oxyJam_RM.out& 

nohup RepeatMasker -species chicken -gff -dir stiNae_RM/ stiNae.1.fasta.gz &> nohup_stiNae_RM.out& 

