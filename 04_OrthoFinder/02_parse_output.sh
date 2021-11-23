# in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/

for file in results/Results_Nov09/Orthogroup_Sequences/*.fa;
do
  export n=$(grep -c '^>' $file)
  echo $file'_'$n >> seq_counts.txt
done

Rscript clean_ogs.R
