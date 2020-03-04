# extract CNEEs from HAL alignment 
cd ../Comparative-Annotation-Toolkit
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	halLiftover --outPSLWithName data/galloanserae.hal galGal ../cnees/galGal6_final_merged_CNEEs_named.bed $SP $SP.psl &
done
exit
mv *.psl ../cnees/alignments
cd ../cnees/alignments

# get ratite script from Tim, edit for specific genomes in HAL file and file names
wget https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl
cp ../Comparative-Annotation-Toolkit/data/galloanserae.hal .
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
cpanm Math::Round
perl parse_cnee_halLiftover.pl 

# some qc - check with Tim about what to look for here? 
grep "multiple_liftover_regions" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > multiple_liftovers_byCNEE.log
grep "no_liftover" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > no_liftOver_byCNEE.log


# generating mafft input: bedtools to get a fasta file for each species with an entry for each CNEE, then split and merge by CNEE name, output should be one file per CNEE with fasta header = species (instead of one file per species with fasta header = CNEE)

for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bedtools getfasta -name -fi ../../ducks_cactus/for_cnees/$SP.defline.fasta -bed ${SP}_cnees_parsed_liftover.bed -name -s > $SP.cnees.fa &
done

#use bioawk to fix up - kind of janky
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' $SP.cnees.fa >> all_cnees.tab
done

#check
cut -f1,1 all_cnees.tab | sort | uniq -c
333952 anaPla
347440 ansBra
346290 ansCyg
346839 ansInd
347588 braCan
330133 colVir
346892 cotJap
356479 galGal
342082 hetAtr
342247 netAur
350558 numMel
341976 oxyJam
337192 stiNae
345733 syrMik
350140 tymCupPin

wc -l *_liftover.bed
333952 anaPla_cnees_parsed_liftover.bed
347440 ansBra_cnees_parsed_liftover.bed
346290 ansCyg_cnees_parsed_liftover.bed
346839 ansInd_cnees_parsed_liftover.bed
347588 braCan_cnees_parsed_liftover.bed
330133 colVir_cnees_parsed_liftover.bed
346892 cotJap_cnees_parsed_liftover.bed
356479 galGal_cnees_parsed_liftover.bed
342082 hetAtr_cnees_parsed_liftover.bed
342247 netAur_cnees_parsed_liftover.bed
350558 numMel_cnees_parsed_liftover.bed
341976 oxyJam_cnees_parsed_liftover.bed
337192 stiNae_cnees_parsed_liftover.bed
345733 syrMik_cnees_parsed_liftover.bed
350140 tymCupPin_cnees_parsed_liftover.bed
5165541 total

# convert back to fasta by column 2 
mkdir unaligned
awk '{print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cnees.tab

mkdir -p aligned
cut -f4,4 galGal_cnees_parsed_liftover.bed | split -a 3 -d  -l 1000 - batch
mv batch* aligned

