######################################
## Running scripts for fastq to VCF ##
######################################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/CompPopGen/duck_popgen

python3 ../../00_setup_local_fastq.py --config REUducks/netAur_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/CompPopGen/duck_popgen/REUducks --abbv Nauritus --config_out REUducks/netAur_config.txt
python3 ../../00_setup_local_fastq.py --config REUducks/oxyJam_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/CompPopGen/duck_popgen/REUducks --abbv Ojamaicensis --config_out REUducks/oxyJam_config.txt
python3 ../../00_setup_local_fastq.py --config REUducks/stiNae_REU_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/CompPopGen/duck_popgen/REUducks --abbv Snaevosa_REU --config_out REUducks/stiNae_REU_config.txt

# manually add in required options to new configs before running 01_download_qc.py
python3 ../../01_download_qc.py --config REUducks/netAur_config.txt
python3 ../../01_download_qc.py --config REUducks/oxyJam_config.txt
python3 ../../01_download_qc.py --config REUducks/stiNae_REU_config.txt

#conda create -n mpl python=3.6 anaconda matplotlib bedtools vcftools 
conda activate mpl

python3 ../../02_dedup_gather_metrics.py --config REUducks/netAur_config.txt
python3 ../../02_dedup_gather_metrics.py --config REUducks/oxyJam_config.txt
python3 ../../02_dedup_gather_metrics.py --config REUducks/stiNae_REU_config.txt

# based on recommendation from 02_dedup_gather_metrics.py, manually change pipeline version in config (or add option when running 03_haplotypecalling.py) 

python3 ../../03_haplotypecalling.py --config REUducks/netAur_config.txt
python3 ../../03_haplotypecalling.py --config REUducks/oxyJam_config.txt
python3 ../../03_haplotypecalling.py --config REUducks/stiNae_REU_config.txt

python3 ../../04_gvcf_filtering.py --config REUducks/netAur_config.txt
python3 ../../04_gvcf_filtering.py --config REUducks/oxyJam_config.txt
python3 ../../04_gvcf_filtering.py --config REUducks/stiNae_REU_config.txt

python3 ../../05_calculate_coverage.py --config REUducks/netAur_config.txt
python3 ../../05_calculate_coverage.py --config REUducks/oxyJam_config.txt
python3 ../../05_calculate_coverage.py --config REUducks/stiNae_REU_config.txt
