## PhyloAcc  
  
This directory contains:  
* scripts `01*_prepPhyloAcc_top*.sh`, which contain the code required to prepare the input data for PhyloAcc for each of the three possible phylogenetic topologies  
* script `02_run_phyloAcc.sh`, which contains the `sbatch` commands to submit the PhyloAcc scripts for each topology  
* scripts `03_parse_phyloAcc.sh`, `04_cleanup_phyloAcc_ncbi.R` and `05_prelim_phyloAcc_ncbi.R`, which contain the code to parse the PhyloAcc output for each topology  
* scripts `06_phyloP_cleanup_cnees.sh` and `07_phyloP_cleanup_cnees.R`, which contains the code to remove all CNEEs that don't show strong evidence for conservation in Galloanserae from phyloP under all topologies after FDR correction  
* `phyloAcc_allDucks` directory, which contains the code and outputs for identical PhyloAcc runs using each focal species as the target relative to the rest of the Galloanserae species used in the whole genome alignment  
* `phyloAcc_control` directory, which contains the code and outputs for an identical PhyloAcc run using all Anseriform species as the target relative to the Galliform species of the whole genome alignment  
  
  
Data related to these analyses are archived in `/n/holylfs05/LABS/informatics/Users/swuitchik/ducks/02_ncbi_analyses/03_cnees`
