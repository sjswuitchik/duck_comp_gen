## Concat output from PhyloAcc ## 

VERSION=$1
cd top1_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top1_${VERSION}.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top1_${VERSION}.txt
cd ../top2_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top2_${VERSION}.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top2_${VERSION}.txt
cd ../top3_outs
../concat_output.sh elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top3_${VERSION}.txt
../concat_output.sh  rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top3_${VERSION}.txt
cd ..
