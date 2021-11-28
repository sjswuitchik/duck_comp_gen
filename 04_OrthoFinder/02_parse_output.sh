# in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/

for file in results/Results_Nov09/Orthogroup_Sequences/*.fa;
do
  export n=$(grep -c '^>' $file)
  echo $file'_'$n >> seq_counts.txt
done

Rscript clean_ogs.R

cd /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/
for file in anaPla ansCyg colVir cotJap galGal numMel
do
  grep -o 'protein_id=.._[0-9]*.[0-9]' $file/${file}_cds_from_genomic.fna > ${file}_protID.tsv
done

Rscript clean_prots.R
