## in /n/holyscratch01/informatics/swuitchik/ducks/compGen

# copy over single copy OG seqs
mkdir -p single_copies/fastas
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Single_Copy_Orthologue_Sequences/*.fa single_copies/fastas
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Orthogroups/Orthogroups_SingleCopyOrthologues.txt single_copies
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Orthogroups/Orthogroups.tsv single_copies
cp -vr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09/Species_Tree/SpeciesTree_rooted.txt single_copies

cd single_copies

# write species-specific clean OG files
Rscript singleCopy_clean_ogs.R

# 
