## MK & SnIPRE tests  
  
This directory contains analyses following an earlier version of the MK pipeline detailed in https://github.com/sjswuitchik/compPopGen_ms/tree/master/MKpipeline. The files and subdirectories contained here are:  
  
* script `00_setup_pipeline.sh`, which contains the necessary pre-processing steps 
* script `01_pipeline.sh`, which contains the semi-automated pipeline steps
* script `02_prelim_results.R`, which contains the code required to produce the outputs
* `outputs` directory, which contains the oupts from `01_pipeline.sh` and `02_prelim_results.R`
* `required_files` directory, which contains the sources files for SnIPRE
* `subscripts` directory, which contains any scripts used in the main scripts detailed above
* `directory_tree.pdf`, which outlines how the working directories should be set up for `02_pipeline.sh` to automatically navigate the directory structure
  
  
Data associated with these analyses are archived at `/n/holylfs05/LABS/informatics/Users/swuitchik/ducks/02_ncbi_analyses/06_compGen/mk_tests`
