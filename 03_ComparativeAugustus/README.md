## Genome annotation using Comparative Augustus 
  
This directory contains:  
* script `01a_run_CompAug_denovo.sh`, which contains the code for genome annotation using Comparative Augustus _de novo_ 
* script `01b_run_CompAug_hints.sh`, which contains the code for genome annotation using Comparative Augustus with RNA seq hints from chicken, mallard, and ruddy duck
* script `02_BUSCO.sh`, which contains the code for quality checking the annotations with BUSCO
* `input_files` directory, which contains files generated in the Comparative Augustus & BUSCO scripts that are required as inputs in subsequent steps of the analyses (e.g., hints.tbl)  or are specific input files distributed with the software that have been edited (e.g., config.ini)
* `subscripts` directory, which contains 
