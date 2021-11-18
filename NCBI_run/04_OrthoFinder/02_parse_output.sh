# in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Orthogroup_Sequences

for file in *.fa;
do
  export n=$(grep -c '^>' $file)
  echo $file'_'$n >> seq_counts.txt
done

Rscript clean_ogs.R
