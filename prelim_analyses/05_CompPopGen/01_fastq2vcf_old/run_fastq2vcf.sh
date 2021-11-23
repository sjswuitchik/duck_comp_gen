######################################
## Running scripts for fastq to VCF ##
######################################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen
python3 00_setup_local_fastq.py --config configs/hetAtr_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen --abbv Hatricapilla --config_out configs/hetAtr_config.txt
python3 00_setup_local_fastq.py --config configs/stiNae_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen --abbv Snaevosa --config_out configs/stiNae_config.txt

# manually add in required options to new configs before running 01_download_qc.py

#conda create -n mpl python=3.6 anaconda matplotlib bedtools vcftools 
conda activate mpl

python3 01_download_qc.py --config configs/hetAtr_config.txt
python3 01_download_qc.py --config configs/stiNae_config.txt

python3 02_dedup_gather_metrics.py --config configs/hetAtr_config.txt 
python3 02_dedup_gather_metrics.py --config configs/stiNae_config.txt

# based on recommendation from 02_dedup_gather_metrics.py, manually change pipeline version in config (or add option when running 03_haplotypecalling.py) 

python3 03_haplotypecalling.py --config configs/hetAtr_config.txt 
python3 03_haplotypecalling.py --config configs/stiNae_config.txt

python3 04_genotypegvcf_filtering.py --config configs/hetAtr_config.txt 
python3 04_genotypegvcf_filtering.py --config configs/stiNae_config.txt

python3 05_calculate_coverage.py --config configs/hetAtr_config.txt 
python3 05_calculate_coverage.py --config configs/stiNae_config.txt

