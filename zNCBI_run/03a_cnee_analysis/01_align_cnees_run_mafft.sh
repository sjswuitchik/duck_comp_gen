###########################
## CNEE alignment - NCBI ##
###########################

module load Anaconda3/2019.10

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis

mkdir -p alignments
cp /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/02_wga/gallo_ncbi.hal .
# CNEE set was already made from earlier analysis, so copy over (from duck_comp_gen/03a_cnee_analysis/01_assemble_cnees.sh)
cp /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/cnees/galGal6_final_merged_CNEEs_named.bed .

singularity shell --cleanenv /n/singularity_images/informatics/cat/cat:20200116.sif
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap galGal hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin;
do
	halLiftover --outPSLWithName gallo_ncbi.hal galGal galGal6_final_merged_CNEEs_named.bed $SP $SP.psl &
done
exit

mv *.psl alignments/
cd alignments/

# get ratite script from Tim, edit for specific genomes in HAL file and file names
wget https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl
conda activate py27
#cpanm Math::Round
perl parse_cnee_halLiftover.pl 

# some qc
grep "multiple_liftover_regions" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > multiple_liftovers_byCNEE.log
grep "no_liftover" galGal6_final_merged_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > no_liftOver_byCNEE.log

# generating mafft input: bedtools to get a fasta file for each species with an entry for each CNEE, then split and merge by CNEE name, output should be one file per CNEE with fasta header = species (instead of one file per species with fasta header = CNEE)

for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bedtools getfasta -name -fi /n/holylfs/LABS/informatics/swuitchik/ducks/ducks_cactus/for_cnees/$SP.defline.fasta -bed ${SP}_cnees_parsed_liftover.bed -name -s > $SP.cnees.fa &
done

#use bioawk to fix up - kind of janky
for SP in anaPla ansBra ansCyg ansInd braCan colVir cotJap hetAtr netAur numMel oxyJam stiNae syrMik tymCupPin
do
	bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' $SP.cnees.fa >> all_cnees.tab
done

#check
cut -f1,1 all_cnees.tab | sort | uniq -c




wc -l *_liftover.bed





# convert back to fasta by column 2 
mkdir unaligned
awk '{print ">"$1 >> "unaligned/"$2".fa"; print $3 >> "unaligned/"$2".fa"; close("unaligned/"$2".fa")}' all_cnees.tab

mkdir -p aligned
cut -f4,4 galGal_cnees_parsed_liftover.bed | split -a 3 -d  -l 1000 - batch
mv batch* aligned/

# run MAFFT 
cd aligned/
sbatch run_mafft.sh

