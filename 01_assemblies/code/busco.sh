## Running BUSCO on unmasked assemblies

module load centos6 BUSCO/3.0.2-fasrc01 python/3.6.0-fasrc01
tar -xzf aves_odb9.tar.gz
cp /n/sw/fasrcsw/apps/Core/BUSCO/3.0.2-fasrc01/bin/../config/config.ini .

# need to edit config file with specific paths for dependencies

whereis makeblastdb
whereis tblasn
whereis hmmsearch
whereis augustus #this path has to be pasted four times in the appropriate places in the config - once with /bin at the end for the main path, and three times with /scripts instead of /bin at the end of the path for each script
whereis etraining
whereis Rscript # only needed if you're going to use the plot function

cp -r /n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/config/ /scratch/swuitchik/raw_data_ducks/genomes/

export BUSCO_CONFIG_FILE=/scratch/swuitchik/raw_data_ducks/genomes/config.ini
export AUGUSTUS_CONFIG_PATH=/scratch/swuitchik/raw_data_ducks/genomes/config

nohup run_BUSCO.py -i hetAtr.1.fasta -o hetAtr_busco -l aves_odb9 -m geno -f &> nohup_hetAtr_busco.out&  

nohup run_BUSCO.py -i netAur.1.fasta -o netAur_busco -l aves_odb9 -m geno &> nohup_netAur_busco.out&

nohup run_BUSCO.py -i oxyJam.1.fasta -o oxyJam_busco -l aves_odb9 -m geno &> nohup_oxyJam_busco.out&

nohup run_BUSCO.py -i stiNae.1.fasta -o stiNae_busco -l aves_odb9 -m geno &> nohup_stiNae_busco.out&

