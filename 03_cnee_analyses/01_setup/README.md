## Assembling & aligning CNEEs for PhyloAcc input  
  
This directory contains:  
* script `01_assemble_cnees.sh`, which contains the code used to build the consensus set of CNEEs from the published literature on vertebrates and Aves  
* script `02_align_mafft.sh`, which contains the code to liftover CNEEs and generate FASTAs for each species in the whole genome alignment, then align CNEEs using MAFFT  
* script `03_catseq.sh`, which contains the code for concatenating CNEEs & generating a partitions BED  
* script `04_4d_sites.sh`, which contains the code for creating 4d neutral models for each possibly phylogenetic topology  
* `subscripts` directory, which contains any scripts used in the main scripts
