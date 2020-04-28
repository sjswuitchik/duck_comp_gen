# extract CNEEs from HAL alignment 
singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	halLiftover --outPSLWithName galloanserae.hal galGal galGal6_final_merged_CNEEs_named.bed $SP $SP.psl &
done
exit
mv *.psl alignments
cd alignments

# get ratite script from Tim, edit for specific genomes in HAL file and file names
wget https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl
cp ../galloanserae.hal .
cp ../galGal6_final_merged_CNEEs_named.bed .
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
#cpanm Math::Round
perl parse_cnee_halLiftover.pl 

# some qc - check with Tim about what to look for here? 
grep "multiple_liftover_regions" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > multiple_liftovers_byCNEE.log
grep "no_liftover" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > no_liftOver_byCNEE.log


# generating mafft input: bedtools to get a fasta file for each species with an entry for each CNEE, then split and merge by CNEE name, output should be one file per CNEE with fasta header = species (instead of one file per species with fasta header = CNEE)

for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bedtools getfasta -name -fi /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees/$SP.defline.fasta -bed ${SP}_cnees_parsed_liftover.bed -name -s > $SP.cnees.fa &
done

#use bioawk to fix up - kind of janky
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' $SP.cnees.fa >> all_cnees.tab
done

#check
cut -f1,1 all_cnees.tab | sort | uniq -c
336525 anaPla
354452 ansBra
351133 ansCyg
349543 ansInd
350456 braCan
333030 colVir
556034 cotJap
635687 galGal
499698 hetAtr
459919 netAur
529159 numMel
526466 oxyJam
491231 stiNae
349247 syrMik
353303 tymCupPin

wc -l *_liftover.bed
351462 anaPla_cnees_parsed_liftover.bed
365817 ansBra_cnees_parsed_liftover.bed
364498 ansCyg_cnees_parsed_liftover.bed
365169 ansInd_cnees_parsed_liftover.bed
366023 braCan_cnees_parsed_liftover.bed
347679 colVir_cnees_parsed_liftover.bed
365255 cotJap_cnees_parsed_liftover.bed
375591 galGal_cnees_parsed_liftover.bed
360137 hetAtr_cnees_parsed_liftover.bed
360316 netAur_cnees_parsed_liftover.bed
369233 numMel_cnees_parsed_liftover.bed
359926 oxyJam_cnees_parsed_liftover.bed
355035 stiNae_cnees_parsed_liftover.bed
364126 syrMik_cnees_parsed_liftover.bed
368743 tymCupPin_cnees_parsed_liftover.bed
5439010 total

# convert back to fasta by column 2 
mkdir unaligned
awk '{print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cnees.tab

mkdir -p aligned
cut -f4,4 galGal_cnees_parsed_liftover.bed | split -a 3 -d  -l 1000 - batch
mv batch* aligned/

