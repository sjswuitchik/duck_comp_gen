# Comparative genomics of brood parasitism in the black-headed duck

Authors:

Sara Smith Wuitchik (Assistant Professor, Mount Royal University; ssmith6@mtroyal.ca)

LaDeana Hillier (Department of Genome Sciences, Washington University; lhillier@uw.edu)

Chris Balakrishnan (Associate Professor, East Carolina University; balakrishnanc@ecu.edu)

Wes Warren (Professor, University of Missouri; warrenwc@missouri.edu)

Mike Sorenson (Professor, Boston University; msoren@bu.edu)

Tim Sackton (Director of Bioinformatics, Informatics Group, Harvard University; tsackton@g.harvard.edu)



Genome assembly and comparative genomics project with the freckled duck (Stictonetta naevosa; stiNae), ruddy duck (Oxyura jamaicensis; oxyJam), African pygmy-goose (Nettapus auritus; netAur), and black-headed duck (Heteronetta atricapilla; hetAtr) from 10x data using Supernova assembly.

Code is currently being organized and optimized. Code and select data* related to the following analyses can be found in these directories: 

`01_assemblies`: Genome assemblies using Supernova of 10X data and quality checks with BUSCO 
`02_wga`: Generation of a whole genome alignment of Galloanserae genomes using CACTUS  
`03_ComparativeAugustus`: Generation of genome annotations (both _de novo_ and hinted) with Comparative Augustus and quality checks with BUSCO  
`03_cnee_analyses`: Compilation of conserved non-coding elements from Aves and vertebrates and multiple PhyloAcc (https://phyloacc.github.io/) analyses  
`04_OrthoFinder`: Generation of orthogroups using OrthoFinder (https://github.com/davidemms/OrthoFinder)  
`04_polytomy_resolution`: Resolution of phylogenetic polytomy between the focal species using coding and non-coding sequences  
`05_CompPopGen`: VCF generation and quality checks using snpArcher (https://github.com/harvardinformatics/snpArcher), McDonald-Kreitman tests & SnIPRE (Eilertson et al. 2012) for selection using the framework outlined in https://github.com/sjswuitchik/compPopGen_ms/tree/master/MKpipeline, tests for selection using HyPhy (https://github.com/veg/hyphy), demographic inference using Stairway (Liu & Fu 2020), and identification of selective sweeps using SweepFinder2 (DeGiorgio et al. 2016)  
`prelim_analyses`: initial analyses of chromosome-only assemblies; not current. 





\* select data available through GitHub where file sizes are compatible with GitHub file size permissions. Genome assemblies are available on NCBI (BioProject PRJNA588796).
