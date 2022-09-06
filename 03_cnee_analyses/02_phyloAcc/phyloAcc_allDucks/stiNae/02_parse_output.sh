# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/phyloAcc_allDucks/stiNae/PhyloAcc

mkdir slurm_outs
mv slurm-* slurm_outs/

#*_elem_lik is the likelihood, needs to be parsed to remove 0s
#*rate_postZ_M2.txt is posterior probs

## top1
#make headers
head top1_outs/batch000_elem_lik.txt -n 1 > elem_lik_top1.header
head top1_outs/batch000_rate_postZ_M2.txt -n 1 > rate_postZ_M2_top1.header

#headers: No.     ID      loglik_NUll     loglik_RES      loglik_all      log_ratio       loglik_Max1     loglik_Max2     loglik_Max3
cat top1_outs/*_elem_lik.txt | awk 'BEGIN {OFS = "\t"} {if ($3 != 0) {print}}'  | grep -v  "^No" > top1_combined_elem_lik.temp 
cat top1_outs/*_rate_postZ_M2.txt | grep -v  "^No" >  top1_combined_postZ_M2.temp  

cat elem_lik_top1.header top1_combined_elem_lik.temp > top1_combined_elem_lik.txt
cat rate_postZ_M2_top1.header top1_combined_postZ_M2.temp >  top1_combined_postZ_M2.txt

## top2
#make headers
head top2_outs/batch000_elem_lik.txt -n 1 > elem_lik_top2.header
head top2_outs/batch000_rate_postZ_M2.txt -n 1 > rate_postZ_M2_top2.header

#headers: No.     ID      loglik_NUll     loglik_RES      loglik_all      log_ratio       loglik_Max1     loglik_Max2     loglik_Max3
cat top2_outs/*_elem_lik.txt | awk 'BEGIN {OFS = "\t"} {if ($3 != 0) {print}}'  | grep -v  "^No" > top2_combined_elem_lik.temp 
cat top2_outs/*_rate_postZ_M2.txt | grep -v  "^No" >  top2_combined_postZ_M2.temp  

cat elem_lik_top2.header top2_combined_elem_lik.temp > top2_combined_elem_lik.txt
cat rate_postZ_M2_top2.header top2_combined_postZ_M2.temp >  top2_combined_postZ_M2.txt

## top3
#make headers
head top3_outs/batch000_elem_lik.txt -n 1 > elem_lik_top3.header
head top3_outs/batch000_rate_postZ_M2.txt -n 1 > rate_postZ_M2_top3.header

#headers: No.     ID      loglik_NUll     loglik_RES      loglik_all      log_ratio       loglik_Max1     loglik_Max2     loglik_Max3
cat top3_outs/*_elem_lik.txt | awk 'BEGIN {OFS = "\t"} {if ($3 != 0) {print}}'  | grep -v  "^No" > top3_combined_elem_lik.temp 
cat top3_outs/*_rate_postZ_M2.txt | grep -v  "^No" >  top3_combined_postZ_M2.temp  

cat elem_lik_top3.header top3_combined_elem_lik.temp > top3_combined_elem_lik.txt
cat rate_postZ_M2_top3.header top3_combined_postZ_M2.temp >  top3_combined_postZ_M2.txt

gzip *.txt
rm *.temp
