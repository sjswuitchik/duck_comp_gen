# run_supernova.sh 

# useage: nohup ./run_supernova.sh sample path &> nohup.out&

module load supernova/2.1.1-fasrc01
FASTQ="$2"
SAMPLE="$1"

supernova run --id $SAMPLE --fastqs "$FASTQ" --maxreads 448000000 --localcores 24 --localmem 400

