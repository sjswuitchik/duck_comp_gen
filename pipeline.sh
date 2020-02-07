############################
## Genome assembly script ##
############################

# Supernova assembly

FASTQ="$2"
SAMPLE="$1"
supernova run --id $SAMPLE --fastqs "$FASTQ" --maxreads 448000000 --localcores 24 --localmem 400

# Supernova: create output
supernova mkoutput --style=pseudohap2 --asmdir=/scratch/swuitchik/raw_data_ducks/stiNae/outs/assembly --outprefix=stiNae > mkout_stiNae.out

# BUSCO: on unmasked assemblies
# download appropriate database from https://busco.ezlab.org/
tar -xzf tar -xzf aves_odb9.tar.gz
cp /n/sw/fasrcsw/apps/Core/BUSCO/3.0.2-fasrc01/bin/../config/config.ini .

# need to edit config file with specific paths for dependencies

whereis makeblastdb
whereis tblastn
whereis hmmsearch
whereis augustus #this path has to be pasted four times in the appropriate places in the config - once with /bin at the end for the main path, and three times with /scripts instead of /bin at the end of the path for each script
whereis etraining
whereis Rscript # only needed if you're going to use the plot function

cp -r /n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/config/ .
export BUSCO_CONFIG_FILE=/scratch/swuitchik/raw_data_ducks/genomes/config.ini
export AUGUSTUS_CONFIG_PATH=/scratch/swuitchik/raw_data_ducks/genomes/config

run_BUSCO.py -i stiNae.1.fasta -o stiNae_busco -l aves_odb9 -m geno > stiNae_busco.out

# RepeatMasker: on ordered and oriented chromosomal fastas
RepeatMasker -species chicken -xsmall -gff -dir stiNae_RM_chrom/ BPBGC_Snae1.0.chromosomal_ordered_oriented.fasta > stiNae_RM.out
