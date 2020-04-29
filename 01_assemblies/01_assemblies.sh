####################################
## Genome assemblies & initial QC ##
####################################

# Supernova assembly

nohup ./run_supernova.sh stiNae 01_stiNae/ &> nohup_stiNae.out&
nohup ./run_supernova.sh oxyJam 02_oxyJam/ &> nohup_oxyJam.out&
nohup ./run_supernova.sh netAur 03_netAur/ &> nohup_netAur.out&
nohup ./run_supernova.sh hetAtr 04_hetAtr/ &> nohup_hetAtr.out&

# Make Supernova output (both pseudohaplotypes)

module load supernova/2.1.1-fasrc01

supernova mkoutput --style=pseudohap2 --asmdir=/n/holylfs/LABS/informatics/swuitchik/ducks/outputs/stiNae/stiNae_outs/outs/assembly --outprefix=stiNae > mkout_stiNae.out

supernova mkoutput --style=pseudohap2 --asmdir=/n/holylfs/LABS/informatics/swuitchik/ducks/outputs/oxyJam/oxyJam_outs/outs/assembly --outprefix=oxyJam > mkout_oxyJam.out

supernova mkoutput --style=pseudohap2 --asmdir=/n/holylfs/LABS/informatics/swuitchik/ducks/outputs/netAur/netAur_outs/outs/assembly --outprefix=netAur > mkout_netAur.out

supernova mkoutput --style=pseudohap2 --asmdir=/n/holylfs/LABS/informatics/swuitchik/ducks/outputs/hetAtr/hetAtr_outs/outs/assembly --outprefix=hetAtr > mkout_hetAtr.out

# BUSCO for QC 

module load centos6 BUSCO/3.0.2-fasrc01 python/3.6.0-fasrc01
tar -xzf aves_odb9.tar.gz
cp /n/sw/fasrcsw/apps/Core/BUSCO/3.0.2-fasrc01/bin/../config/config.ini .

# need to edit config file with paths to specific dependencies
whereis makeblastdb
whereis tblastn
whereis hmmsearch
whereis augustus #this path has to be pasted four times in the appropriate places in the config - once with /bin at the end for the main path, and three times with /scripts instead of /bin at the end of the path for each script
whereis etraining
whereis Rscript # only needed if you're going to use the plot function

# copy these paths individually into the config file - be sure to not include the command
# eg/ for makeblastdb, the path is /n/sw/fasrcsw/apps/Core/ncbi-blast/2.2.31+-fasrc01/bin/makeblastdb but you only copy /n/sw/fasrcsw/apps/Core/ncbi-blast/2.2.31+-fasrc01/bin/ into the config file

cd /n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/
cp -r config/ /scratch/swuitchik/raw_data_ducks/genomes/
cd /scratch/swuitchik/raw_data_ducks/genome

export BUSCO_CONFIG_FILE=/scratch/swuitchik/raw_data_ducks/genomes/config.ini
export AUGUSTUS_CONFIG_PATH=/scratch/swuitchik/raw_data_ducks/genomes/config

nohup run_BUSCO.py -i hetAtr.1.fasta -o hetAtr_busco -l aves_odb9 -m geno -f &> nohup_hetAtr_busco.out&  

nohup run_BUSCO.py -i netAur.1.fasta -o netAur_busco -l aves_odb9 -m geno &> nohup_netAur_busco.out&

nohup run_BUSCO.py -i oxyJam.1.fasta -o oxyJam_busco -l aves_odb9 -m geno &> nohup_oxyJam_busco.out&

nohup run_BUSCO.py -i stiNae.1.fasta -o stiNae_busco -l aves_odb9 -m geno &> nohup_stiNae_busco.out&

