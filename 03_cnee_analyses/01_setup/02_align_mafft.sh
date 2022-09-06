###########################
## CNEE alignment - NCBI ##
###########################

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis

mkdir -p alignments
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/02_wga/gallo_ncbi.hal .
# CNEE set was already made from earlier analysis, so copy over (from duck_comp_gen/03a_cnee_analysis/01_assemble_cnees.sh)
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/galGal6_final_merged_CNEEs_named.bed .

singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200604.sif
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	halLiftover --outPSLWithName gallo_ncbi.hal galGal galGal6_final_merged_CNEEs_named.bed $SP $SP.psl &
done
exit

mv *.psl alignments/
cp galGal6_final_merged_CNEEs_named.bed alignments/
cd alignments/

# get ratite script from Tim, edit for file names (line 14) and genomes in the HAL (line 39) 
wget https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl
conda activate py27
#cpanm Math::Round
perl parse_cnee_halLiftover.pl 

# some qc
grep "multiple_liftover_regions" final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > multiple_liftovers_byCNEE.log
grep "no_liftover" final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > no_liftOver_byCNEE.log

# generating mafft input: bedtools to get a fasta file for each species with an entry for each CNEE, then split and merge by CNEE name, output should be one file per CNEE with fasta header = species (instead of one file per species with fasta header = CNEE)
cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/02_wga/*.fasta /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/fastas
cd ../fastas/
/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/brename -p ".defline.fasta" -r ".fasta"
/n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cesar/2bitdir/brename -p ".ncbi.fasta" -r ".fasta"
cd ../alignments

for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	bedtools getfasta -name -fi ../fastas/$SP.fasta -bed ${SP}_cnees_parsed_liftover.bed -name -s > $SP.cnees.fa &
done

#use bioawk to fix up - kind of janky but doesn't create a ton of intermediate files
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' $SP.cnees.fa >> all_cnees.tab
done

#check
cut -f1,1 all_cnees.tab | sort | uniq -c
351418 anaPla
365896 ansBra
364527 ansCyg
365247 ansInd
366083 braCan
347707 colVir
365254 cotJap
375591 galGal
361927 hetAtr
361097 netAur
369276 numMel
362552 oxyJam
359132 stiNae
364132 syrMik
368740 tymCupPin

wc -l *_liftover.bed
351418 anaPla_cnees_parsed_liftover.bed
365896 ansBra_cnees_parsed_liftover.bed
364527 ansCyg_cnees_parsed_liftover.bed
365247 ansInd_cnees_parsed_liftover.bed
366083 braCan_cnees_parsed_liftover.bed
347707 colVir_cnees_parsed_liftover.bed
365254 cotJap_cnees_parsed_liftover.bed
375591 galGal_cnees_parsed_liftover.bed
361927 hetAtr_cnees_parsed_liftover.bed
361097 netAur_cnees_parsed_liftover.bed
369276 numMel_cnees_parsed_liftover.bed
362552 oxyJam_cnees_parsed_liftover.bed
359132 stiNae_cnees_parsed_liftover.bed
364132 syrMik_cnees_parsed_liftover.bed
368740 tymCupPin_cnees_parsed_liftover.bed
5448579 total

# convert back to fasta by column 2 
mkdir unaligned
awk '{print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cnees.tab

mkdir -p aligned
cut -f4,4 galGal_cnees_parsed_liftover.bed | split -a 3 -d  -l 1000 - batch
mv batch* aligned/

# run MAFFT 
cd aligned/
sbatch run_mafft.sh

