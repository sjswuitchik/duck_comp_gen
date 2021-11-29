##  in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/

# get OGs that are in hetAtr, in at least 50% of spp, and capped at 25 seq
Rscript clean_ogs.R

## NCBI spp
# get protein IDs for each species on NCBI
cd ncbi_data/
for file in anaPla ansCyg colVir cotJap galGal numMel;
do
  grep -o 'protein_id=.._[0-9]*.[0-9]' $file/${file}_cds_from_genomic.fna > ${file}_protID.tsv
done

# check numbers match
for file in anaPla ansCyg colVir cotJap galGal numMel;
do
  grep -c 'protein_id=' $file/${file}_cds_from_genomic.fna
  wc -l ${file}_protID.tsv
done

#39905
#39905 anaPla_protID.tsv
#31811
#31811 ansCyg_protID.tsv
#17165
#0 colVir_protID.tsv
#40894
#40894 cotJap_protID.tsv
#49681
#49681 galGal_protID.tsv
#43240
#43240 numMel_protID.tsv

# didn't work for colVir, check format of prot ID for that spp
grep -o 'protein_id=OXB[0-9]*.[0-9]' colVir/colVir_cds_from_genomic.fna > colVir_protID.tsv
wc -l colVir_protID.tsv 
# 17165 colVir_protID.tsv

cd ..
for spp in anaPla ansCyg colVir cotJap galGal numMel;
do
  Rscript clean_prots.R ncbi_data/${spp}_protID.tsv ${spp}_cleanID.tsv
done

mkdir clean_ids
mv *_cleanID.tsv clean_ids


## Comp Aug spp
mkdir compAug_data
for spp in ansBra ansInd braCan hetAtr netAur oxyJam stiNae syrMik tymCupPin;
do
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/${spp}.gff compAug_data/
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/genomes/${spp}.ncbi.fasta compAug_data/
  samtools faidx compAug_data/${spp}.ncbi.fasta 
  awk '$3 == "CDS"' compAug_data/${spp}.gff > compAug_data/${spp}.cds.gff
  Rscript 
done




