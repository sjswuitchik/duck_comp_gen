## Concat output from PhyloAcc ## 

# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/PhyloAcc

cd top1_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top1.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top1.txt
cd ../top2_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top2.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top2.txt
cd ../top3_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top3.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top3.txt
