##################################################
## Running scripts for fastq to VCF - REU ducks ##
##################################################

module load Anaconda3/2019.10

# from /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen

python3 00_setup_local_fastq.py --config configs/netAurREU_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen --abbv NauritusREU --config_out configs/netAurREU_config.txt
python3 00_setup_local_fastq.py --config configs/oxyJamREU_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen --abbv OjamaicensisREU --config_out configs/oxyJamREU_config.txt
python3 00_setup_local_fastq.py --config configs/stiNaeREU_local.txt --out_dir /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen --abbv Snaevosa_REU --config_out configs/stiNaeREU_config.txt

# manually add in required options to new configs before running 01_download_qc.py
python3 01_download_qc.py --config configs/netAurREU_config.txt
python3 01_download_qc.py --config configs/oxyJamREU_config.txt
python3 01_download_qc.py --config configs/stiNaeREU_config.txt

#conda create -n mpl python=3.6 anaconda matplotlib bedtools vcftools 
conda activate mpl

python3 02_dedup_gather_metrics.py --config configs/netAurREU_config.txt
python3 02_dedup_gather_metrics.py --config configs/oxyJamREU_config.txt
python3 02_dedup_gather_metrics.py --config configs/stiNaeREU_config.txt

# based on recommendation from 02_dedup_gather_metrics.py, manually change pipeline version in config (or add option when running 03_haplotypecalling.py) 

python3 03_haplotypecalling.py --config configs/netAurREU_config.txt
python3 03_haplotypecalling.py --config configs/oxyJamREU_config.txt
python3 03_haplotypecalling.py --config configs/stiNaeREU_config.txt

python3 04_gvcf_filtering.py --config configs/netAurREU_config.txt
python3 04_gvcf_filtering.py --config configs/oxyJamREU_config.txt
python3 04_gvcf_filtering.py --config configs/stiNaeREU_config.txt

python3 05_calculate_coverage.py --config configs/netAurREU_config.txt
python3 05_calculate_coverage.py --config configs/oxyJamREU_config.txt
python3 05_calculate_coverage.py --config configs/stiNaeREU_config.txt
