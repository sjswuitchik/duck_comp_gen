#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import datetime

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename):
    print("Opening %s"%config_filename)
    config_file = open(config_filename,"r")
    config_info = {}
    sample_ncbi_dict = {}
    sample_ena_dict = {}
    sample_local_dict = {}
    sample_dict = {}
    
    for line in config_file:
        if line[0] == "#":
            pass
        elif line == "\n":
            pass
        else:
            line=line.strip()
            line = line.split(" ")
            if line[0] == "--ABBV":
                config_info["abbv"] = line[1]
            elif line[0] == "--SAMPLE_NCBI":
                if line[1] not in sample_ncbi_dict:
                    sample_ncbi_dict[line[1]] = [line[2]]
                elif line[1] in sample_ncbi_dict:
                    sample_ncbi_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_ncbi_dict"] = sample_ncbi_dict
                config_info["sample_dict"] = sample_dict
            elif line[0] == "--SAMPLE_ENA":
                if line[1] not in sample_ena_dict:
                    sample_ena_dict[line[1]] = [line[2]]
                elif line[1] in sample_ena_dict:
                    sample_ena_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_ena_dict"] = sample_ena_dict
                config_info["sample_dict"] = sample_dict
            elif line[0] == "--SAMPLE_LOCAL":
                if line[1] not in sample_local_dict:
                    sample_local_dict[line[1]] = [line[2]]
                elif line[1] in sample_local_dict:
                    sample_local_dict[line[1]].append(line[2])
                if line[1] not in sample_dict:
                    sample_dict[line[1]] = [line[2]]
                elif line[1] in sample_dict:
                    sample_dict[line[1]].append(line[2])
                config_info["sample_local_dict"] = sample_local_dict
                config_info["sample_dict"] = sample_dict

            elif line[0] == "--GENOME_NCBI":
                config_info["genome_ncbi"] = line[1]
            elif line[0] == "--GENOME_LOCAL":
                config_info["genome_local"] = line[1]
            elif line[0] == "--SAMPLE_LOCAL":
                sys.exit("Local sample files are not supported at this time")
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
    config_file.close()
    
    #Add the sample_ncbi_dict and sample_ena_dict if they don't exist
    if "sample_ncbi_dict" not in config_info:
        config_info["sample_ncbi_dict"] = {}
    
    if "sample_ena_dict" not in config_info:
        config_info["sample_ena_dict"] = {}
    
    if "sample_local_dict" not in config_info:
        config_info["sample_local_dict"] = {}
    
    #Make sure all necessary inputs are present
    try:
        config_info["abbv"]
    except NameError:
        sys.exit("Oops, you forgot to specify a species abbreviation with --ABBV")
    
    try:
        config_info["out_dir"]
    except NameError:
        config_info["out_dir"] = "."
        
    if "genome_ncbi" not in config_info and "genome_local" not in config_info:
        sys.exit("Oops, you forgot to specify a reference genome!")
        
    if len(config_info["sample_dict"]) == 0:
        sys.exit("Oops, you forgot to specify samples!")
    
    #Return objects    
    return(config_info)


#Check if a directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):

        os.makedirs(dir)

#Create generic slurm script
def script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -t {time}\n#SBATCH --mem {mem}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%j.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)

    
#Submit filename to slurm with sbatch, returns job id number
def sbatch_submit(filename):
    proc = Popen('sbatch %s'%filename,shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    stdout = stdout.decode("utf-8","ignore")
    stdout = stdout.strip()
    stdout = stdout.strip('Submitted batch job ')
    return(stdout)


#Check job status of specific jobid: returns job status
def jobid_status(jobid,date):
    proc = Popen('sacct --format state --noheader -j %d -S %s'%(int(jobid),date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split()
    return(lines[0].decode("utf-8","ignore"))
    
    
#Check the status of all jobs, returns dictionary of jobid:status. Ignores jobs with ".ba+" and ".ex+" file extensions.
def all_jobs_status(date):
    proc = Popen('sacct --format jobid,state --noheader -S %s'%(date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    stdout=stdout.decode("utf-8","ignore")
    stderr=stderr.decode("utf-8","ignore")
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split("\n")
    status_dict = {}
    for line in lines:
        line = line.split()
        if len(line) > 1:
            if not re.search('.ba+',line[0]) and not re.search('.ex+',line[0]):
                status_dict[line[0]]=line[1]
    return(status_dict)
    

def num_pend_run(job_id_list,date):
    count = 0
    status_dict = all_jobs_status(date)
    for job in job_id_list:
        if status_dict[job] == "PENDING" or status_dict[job] == "RUNNING":
            count += 1
    return(count)

####Create SBATCH files
#Use FTP to download species genome fasta file. Will obtain FTP location by first downloading the current genbank assembly summary report. Returns name of sbatch script created.
def get_ncbi_genome(sp_dir,fasta_ftp,sp_abbr):

    #Recreate genome directory
    genome_dir = "%s/genome"%(sp_dir)
    
    print("Downloading genome from %s and building relevant indexes"%fasta_ftp)
    
    #Get relevant directories and create filenames for other desired files
    fasta_ftp_split = fasta_ftp.split("/")

    ftp_dir = "/".join(fasta_ftp_split[:-1])
    genome_filename = fasta_ftp_split[-1]
    file_base = fasta_ftp_split[-2]
    
    slurm_script = script_create()
    
    #Load modules, also print samtools and bwa versions
    cmd_1 = 'module load samtools/1.5-fasrc02\nmodule load bwa/0.7.15-fasrc02'
    
    cmd_2 = 'wget -P %s %s'%(genome_dir,fasta_ftp)
    cmd_3 = 'gunzip %s/%s'%(genome_dir,genome_filename)
    cmd_4 = 'mv %s/%s %s/%s.fa'%(genome_dir,genome_filename[:-3],genome_dir,sp_abbr)
    cmd_5 = 'samtools faidx %s/%s.fa'%(genome_dir,sp_abbr)
    cmd_6 = 'bwa index %s/%s.fa'%(genome_dir,sp_abbr)
    cmd_7 = 'samtools dict -o %s/%s.dict %s/%s.fa'%(genome_dir,sp_abbr,genome_dir,sp_abbr)
    cmd_8 = 'wget -P %s %s/%s_assembly_stats.txt'%(genome_dir,ftp_dir,file_base)
    cmd_9 = 'wget -P %s %s/%s_assembly_report.txt'%(genome_dir,ftp_dir,file_base)
    cmd_10 = 'wget -P %s %s/%s_genomic.gff.gz'%(genome_dir,ftp_dir,file_base)
    cmd_11 = 'wget -P %s %s/md5checksums.txt'%(genome_dir,ftp_dir)
    
    cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7,cmd_8,cmd_9,cmd_10,cmd_11]
    
    final_cmd = "\n\n".join(cmd_list)   
    
    #Format sbatch script
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_DL_Index",sp_dir=sp_dir,cmd=final_cmd)

    out_filename = "%s/scripts/01_genome_download_index_%s.sbatch"%(sp_dir,sp_abbr)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close
    
    return(out_filename)


#Make sure that genome and genome indexes are set up in genome directory for mapping. Returns name of sbatch script created.
def process_local_genome(sp_dir,genome_local,sp_abbr,genome_present,bwa_index_present,faidx_index_present,dict_index_present):
    genome_dir = "%s/genome"%(sp_dir)
    
    #Load modules, also print samtools and bwa versions
    cmd_1 = 'module load samtools/1.5-fasrc02\nmodule load bwa/0.7.15-fasrc02'
    
    #Copy genome if necessary
    if genome_present == False:
        print("\nCopying %s to %s/%s.fa"%(genome_local,genome_dir,sp_abbr))
        cmd_2 = 'cp %s %s/%s.fa'%(genome_local,genome_dir,sp_abbr)
    else:
        cmd_2 = ''
    
    #Index with samtools faidx if necessary
    if faidx_index_present == False:
        print("\nIndexing %s/%s.fa with Samtools faidx"%(genome_dir,sp_abbr))
        cmd_3 = 'samtools faidx %s/%s.fa'%(genome_dir,sp_abbr)
    else:
        cmd_3 = ''
        
    #Index with BWA if necessary
    if bwa_index_present == False:
        print("\nIndexing %s/%s.fa with bwa index"%(genome_dir,sp_abbr))
        cmd_4 = 'bwa index %s/%s.fa'%(genome_dir,sp_abbr)
    else:
        cmd_4 = ''
        
    if dict_index_present == False:
        print("\nCreating sequence dictionary file for %s/%s.fa"%(genome_dir,sp_abbr))
        cmd_5 = 'samtools dict -o %s/%s.dict %s/%s.fa'%(genome_dir,sp_abbr,genome_dir,sp_abbr)
    else:
        cmd_5 = ''
        
    final_cmd = "%s\n\n%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4,cmd_5)
    
    #Format sbatch script and write file
    slurm_script = script_create()
    genome_script = slurm_script.format(partition="shared",time="0-8:00",mem="8000",cores="1",nodes="1",jobid="Genome_CP_Index",sp_dir=sp_dir,cmd=final_cmd)

    out_filename = "%s/scripts/02_genome_cp_index_%s.sbatch"%(sp_dir,sp_abbr)
    out_file = open(out_filename,"w")
    out_file.write(genome_script)
    out_file.close

    return(out_filename)

#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def ncbi_sra_download_sbatch(sp_dir,sample_ncbi_dict):
    slurm_script = script_create()
    sra_dl_sbatch_filenames = []
    
    path_to_sratools = "/n/home13/ashultz/sw/progs/sratoolkit.2.8.2-1-centos_linux64/bin/"
    
    '''
    #For now ditching prefetch and using wget to donwload due to some weird ascp errors
    #Paths to various required software
    path_to_ascp="/n/home13/ashultz/.aspera/connect/bin/ascp"
    path_to_ascp_openssh="/n/home13/ashultz/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    path_to_sra_dls = "/n/holylfs/LABS/informatics/ashultz/CompPopGen/raw_sra_files"
    
    for sample in sample_ncbi_dict.keys():
        for sra in sample_ncbi_dict[sample]:
            print('Will download %s for sample %s, split, and run through fastqc'%(sra,sample))
    
            #Load modules and get versions for all programs used
            cmd_1 = 'module load fastqc/0.11.5-fasrc01'
            cmd_2 = r'%sprefetch --force all --max-size 500000000 -a "%s|%s" --ascp-options "-QT -l 10G" %s'%(path_to_sratools,path_to_ascp,path_to_ascp_openssh,sra)
            cmd_3 = r'%sfastq-dump --outdir %s/fastq --gzip --split-files %s/sra/%s.sra'%(path_to_sratools,sp_dir,path_to_sra_dls,sra)
            cmd_4 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra) 
            
            final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/02_sra_download_parse_%s.sbatch"%(sp_dir,sra)
            out_file = open(out_filename,"w")
            out_file.write(sra_script)
            out_file.close
            sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)
    '''       
    
    for sample in sample_ncbi_dict.keys():
        for sra in sample_ncbi_dict[sample]:
            #First check if fastq file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
            fastq_1_filename = '%s/fastq/%s_1.fastq.gz'%(sp_dir,sra)
            bam_filename = '%s/alignment/%s.sorted.bam'%(sp_dir,sra)
            dedup_filename = '%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)
            if os.path.isfile(fastq_1_filename) or os.path.isfile(bam_filename) or os.path.isfile(dedup_filename):
                print('%s_1.fastq.gz or %s.sorted.bam or %s.dedup.sorted.bam already present, skipping'%(sra,sra,sample))
            else:
                print('Will download %s for sample %s, split, and run through fastqc'%(sra,sample))
    
                #Load modules and get versions for all programs used
                cmd_1 = 'module load fastqc/0.11.5-fasrc01'
                cmd_2 = 'wget -P %s/sra/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra'%(sp_dir,sra[0:3],sra[0:6],sra,sra)
                cmd_3 = r'%sfastq-dump --outdir %s/fastq --gzip --split-files %s/sra/%s.sra'%(path_to_sratools,sp_dir,sp_dir,sra)
                
                #add check to see if SRA download worked, and if not, try ENA
                cmd_4 = 'if [ -f %s/fastq/%s_1.fastq.gz ];\n\nthen'%(sp_dir,sra)
                                   
                cmd_5 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra)
                
                cmd_6 = 'else'
                
                if len(sra) < 10:
                    cmd_7 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_1.fastq.gz'%(sp_dir,sra[0:6],sra,sra)
                else:
                    cmd_7 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/%s_1.fastq.gz'%(sp_dir,sra[0:6],sra[-1],sra,sra)
                if len(sra) < 10:
                    cmd_8 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_2.fastq.gz'%(sp_dir,sra[0:6],sra,sra)
                else:
                    cmd_8 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/%s_2.fastq.gz'%(sp_dir,sra[0:6],sra[-1],sra,sra)
                
                cmd_9 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra)
                
                cmd_10 = 'fi'
                
                cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7,cmd_8,cmd_9,cmd_10]
            
                final_cmd = "\n\n".join(cmd_list)
    
        #Format sbatch script
                sra_script = slurm_script.format(partition="shared",time="2-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=sp_dir,cmd=final_cmd)
                out_filename = "%s/scripts/02_sra_download_parse_%s.sbatch"%(sp_dir,sra)
                out_file = open(out_filename,"w")
                out_file.write(sra_script)
                out_file.close
                sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)

def ena_sra_download_sbatch(sp_dir,sample_ena_dict):
    slurm_script = script_create()
    sra_dl_sbatch_filenames = []
    
    path_to_sratools = "/n/home13/ashultz/sw/progs/sratoolkit.2.8.2-1-centos_linux64/bin/"
    
        
    for sample in sample_ena_dict.keys():
        for sra in sample_ena_dict[sample]:
            #First check if fastq file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
            fastq_1_filename = '%s/fastq/%s_1.fastq.gz'%(sp_dir,sra)
            bam_filename = '%s/alignment/%s.sorted.bam'%(sp_dir,sra)
            if os.path.isfile(fastq_1_filename) or os.path.isfile(bam_filename):
                print('%s_1.fastq.gz or %s.sorted.bam already present, skipping'%(sra,sra))
            else:
                print('Will download %s for sample %s, split, and run through fastqc'%(sra,sample))
    
                #Load modules and get versions for all programs used
                cmd_1 = 'module load fastqc/0.11.5-fasrc01'
                
                #ENA uses different directory tree if SRA #s >=0 integers or < 7 integers, so have to take both of those into account.
                if len(sra) < 10:
                    cmd_2 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_1.fastq.gz'%(sp_dir,sra[0:6],sra,sra)
                else:
                    cmd_2 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/%s_1.fastq.gz'%(sp_dir,sra[0:6],sra[-1],sra,sra)
                if len(sra) < 10:
                    cmd_3 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s_2.fastq.gz'%(sp_dir,sra[0:6],sra,sra)
                else:
                    cmd_3 = 'wget -P %s/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/00%s/%s/%s_2.fastq.gz'%(sp_dir,sra[0:6],sra[-1],sra,sra)
                cmd_4 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra) 
            
                final_cmd = "%s\n\n%s\n\n%s\n\n%s"%(cmd_1,cmd_2,cmd_3,cmd_4)
    
    
        #Format sbatch script
                sra_script = slurm_script.format(partition="shared",time="3-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=sp_dir,cmd=final_cmd)
                out_filename = "%s/scripts/02_sra_download_parse_%s.sbatch"%(sp_dir,sra)
                out_file = open(out_filename,"w")
                out_file.write(sra_script)
                out_file.close
                sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)

def process_local_fastq_sbatch(sp_dir,sample_local_dict):
    slurm_script = script_create()
    sra_dl_sbatch_filenames = []
    
        
    for sample in sample_local_dict.keys():
        for sra in sample_local_dict[sample]:
            #First check if fastq file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
            fastq_1_filename = '%s/fastq/%s_1.fastq.gz'%(sp_dir,sra)
            bam_filename = '%s/alignment/%s.sorted.bam'%(sp_dir,sra)
            dedup_filename = '%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)
            if os.path.isfile(fastq_1_filename) or os.path.isfile(bam_filename) or os.path.isfile(dedup_filename):
                print('%s_1.fastq.gz or %s.sorted.bam or %s.dedup.sorted.bam already present, skipping'%(sra,sra,sample))
            else:
                print('Will run %s for sample %s through fastqc'%(sra,sample))
    
                #Load modules and get versions for all programs used
                cmd_1 = 'module load fastqc/0.11.5-fasrc01'
                
                #ENA uses different directory tree if SRA #s >=0 integers or < 7 integers, so have to take both of those into account.
                cmd_2 = 'fastqc -o %s/fastqc %s/fastq/%s_1.fastq.gz %s/fastq/%s_2.fastq.gz'%(sp_dir,sp_dir,sra,sp_dir,sra) 
            
                final_cmd = "%s\n\n%s"%(cmd_1,cmd_2)
    
    
        #Format sbatch script
                sra_script = slurm_script.format(partition="shared",time="3-0:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=sp_dir,cmd=final_cmd)
                out_filename = "%s/scripts/02_process_local_fastq_%s.sbatch"%(sp_dir,sra)
                out_file = open(out_filename,"w")
                out_file.write(sra_script)
                out_file.close
                sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)


#Create sbatch files for trimming, mapping, sorting, indexing and alignment stat creation for each set of paired SRA fastq files.
def fastq_trim_align_stats(sp_dir,sra,sp_abbr,sample):
	print('Will trim, map, sort, and index for %s from sample %s'%(sra,sample))

	#Set read group info: ID = SRA number, SM = sample id, PU: SRA.sample, LB: sample id, PL: illumina (assumes all illumina data.) 
	read_group_info = '@RG\\tID:%s\\tSM:%s\\tPU:%s.%s\\tLB:%s\\tPL:illumina'%(sra,sample,sra,sample,sample)

	#Load necessary modules
	cmd_1 = 'module load NGmerge/0.2-fasrc01\nmodule load bwa/0.7.15-fasrc02\nmodule load jdk/1.8.0_45-fasrc01'

	#If paired end data:
	#Trim adapters with NGmerge
	cmd_2 = 'NGmerge -1 %s/fastq/%s_1.fastq.gz -2 %s/fastq/%s_2.fastq.gz -a -v -o %s/fastq/%s_trimmed -n 8 -z -u 60'%(sp_dir,sra,sp_dir,sra,sp_dir,sra)

	#Map to genome with BWA mem
	cmd_3 = r"bwa mem -M -t 8 -R '%s' %s/genome/%s.fa %s/fastq/%s_trimmed_1.fastq.gz %s/fastq/%s_trimmed_2.fastq.gz > %s/alignment/%s_bwa.sam"%(read_group_info,sp_dir,sp_abbr,sp_dir,sra,sp_dir,sra,sp_dir,sra)

	cmd_4 ='gatk --java-options "-Xmx8g -XX:ParallelGCThreads=6" SortSam -I %s/alignment/%s_bwa.sam -O %s/alignment/%s.sorted.bam --SORT_ORDER=coordinate --CREATE_INDEX=true --COMPRESSION_LEVEL 5'%(sp_dir,sra,sp_dir,sra)

	cmd_5 = 'if [ -f %s/alignment/%s.sorted.bai ]\nthen\n rm %s/alignment/%s_bwa.sam \nfi'%(sp_dir,sra,sp_dir,sra)

	cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5]

	final_cmd = "\n\n".join(cmd_list)

	#Format sbatch script and write file
	slurm_script = script_create()
	genome_script = slurm_script.format(partition="shared",time="4-0:00",mem="36000",cores="8",nodes="1",jobid="Trim_Map",sp_dir=sp_dir,cmd=final_cmd)

	out_filename = "%s/scripts/03_trim_map_stats_%s.sbatch"%(sp_dir,sra)
	out_file = open(out_filename,"w")
	out_file.write(genome_script)
	out_file.close

	return(out_filename)





def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    now = datetime.datetime.now()
    print('Staring work: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get Sample, SRA and Genome attributes
#     '''
#     Example config file format (separate by spaces), only include one genome option: 
#     --ABBV <Species_Abbr>
#     --OUT_DIR <output directory>
#     --SAMPLE_NCBI <SAMPLE_1> <SRA_ID1>
#     --SAMPLE_NCBI <SAMPLE_1> <SRA_ID2>
#     --SAMPLE_NCBI <SAMPLE_2> <SRA_ID1>
#     --SAMPLE_ENA <SAMPLE_1> <SRA_ID1>
#     --SAMPLE_LOCAL <SAMPLE_1> <SRA_ID> (any equivalent ID where the fastq is named SRA_ID_1.fastq.gz and SRA_ID_2.fastq.gz
#     --GENOME_NCBI <NCBI Genome Accession>
#     --GENOME_LOCAL <Local genome fasta file>
#     '''
    
    config_info = extract_config(config_filename)

    #####Check if species directory, logs, scripts, fastq, fastqc, genome, alignment, and stats directories, if not creates them
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("\nOutput will be written to %s\n"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    sra_dir = "%s/sra"%(sp_dir)
    fastq_dir = "%s/fastq"%(sp_dir)
    fastqc_dir = "%s/fastqc"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    alignment_dir = "%s/alignment"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    
    directory_create(sp_dir)
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(sra_dir)
    directory_create(fastq_dir)
    directory_create(fastqc_dir)
    directory_create(genome_dir)
    directory_create(alignment_dir)
    directory_create(stats_dir)
    
    
    #####Prepare genome
    #Create sbatch file to download genome if not already present (checks for abbv.fa), mv to abbv.fa.
    
    genome_path = "%s/%s.fa"%(genome_dir,config_info["abbv"])
    index_path_bwa = "%s/%s.fa.bwt"%(genome_dir,config_info["abbv"])
    index_path_faidx = "%s/%s.fa.fai"%(genome_dir,config_info["abbv"])
    index_path_dict = "%s/%s.dict"%(genome_dir,config_info["abbv"])
    genome_jobnum = []
    
    genome_job_id = None
    
    if "genome_ncbi" in config_info:
        #Check if genome, BWA indexes, and faidx indexes already exist in genome directory. If not, create script to copy and index.
        genome_present = os.path.isfile(genome_path)
        bwa_index_present = os.path.isfile(index_path_bwa)
        faidx_index_present = os.path.isfile(index_path_faidx)
        dict_index_present = os.path.isfile(index_path_dict)

        #Create sbatch script if any missing elements (genome, faidx or bwa index)
        if genome_present == False or bwa_index_present == False or faidx_index_present == False or dict_index_present == False:  
            #Create sbatch script
            genome_sbatch_name = get_ncbi_genome(sp_dir,config_info["genome_ncbi"],config_info["abbv"])
            #Submit sbatch script
            genome_job_id = sbatch_submit(genome_sbatch_name)
        
    elif "genome_local" in config_info:
        #Check if genome, BWA indexes, and faidx indexes already exist in genome directory. If not, create script to copy and index.
        genome_present = os.path.isfile(genome_path)
        bwa_index_present = os.path.isfile(index_path_bwa)
        faidx_index_present = os.path.isfile(index_path_faidx)
        dict_index_present = os.path.isfile(index_path_dict)

        #Create sbatch script if any missing elements (genome, faidx or bwa index)
        if genome_present == False or bwa_index_present == False or faidx_index_present == False or dict_index_present == False:  
            genome_sbatch_name = process_local_genome(sp_dir,config_info["genome_local"],config_info["abbv"],genome_present,bwa_index_present,faidx_index_present,dict_index_present)
            
            #Submit sbatch script
            genome_job_id = sbatch_submit(genome_sbatch_name)
    
    #Only check on genome job if actually submitted. 
    if genome_job_id is not None: 
        #Sleep 30 seconds to give job status time to get to sacct before starting to check    
        sleep(30)
    
        dones = ['COMPLETED','CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL']
        #Check job id status of genome job. If not in one of the 'done' job status categories, wait 30 seconds and check again.
        while jobid_status(genome_job_id,start_date) not in dones:
            #print("Genome not done yet")
            sleep(30)
    
        #Check to make sure job completed, and that all necessary files are present. If not, exit and give information.
        genome_job_completion_status = jobid_status(genome_job_id,start_date)
        if genome_job_completion_status != 'COMPLETED':
            sys.exit("There was a problem downloading and indexing the genome. The job exited with status %s. Please diagnose and fix before moving on"%genome_job_completion_status)
        
        genome_present = os.path.isfile(genome_path)
        bwa_index_present = os.path.isfile(index_path_bwa)
        faidx_index_present = os.path.isfile(index_path_faidx)
        dict_index_present = os.path.isfile(index_path_dict)
        if genome_present == False or bwa_index_present == False or faidx_index_present == False or dict_index_present == False:
            sys.exit("The genome job finished but not all files and indexes are present. Please check, create missing indexes, and resubmit with local genome")

        print("\nGenome download and indexing successfully completed\n")
    
    #####Create sbatch files to download SRA files and use fastq-dump to split
    
    #Create sbatch files
    ncbi_sra_dl_sbatch_filenames = []
    ena_sra_dl_sbatch_filenames = []
    local_fastq_process_sbatch_filenames = []
    
    #Create sbatch files for ncbi
    if len(config_info["sample_ncbi_dict"]) > 0:
        ncbi_sra_dl_sbatch_filenames = ncbi_sra_download_sbatch(sp_dir,config_info["sample_ncbi_dict"])
    #Create sbatch files for ena
    if len(config_info["sample_ena_dict"]) > 0:
        ena_sra_dl_sbatch_filenames = ena_sra_download_sbatch(sp_dir,config_info["sample_ena_dict"])
    #Create sbatch files for local fastqs
    if len(config_info["sample_local_dict"]) > 0:
        local_fastq_process_sbatch_filenames = process_local_fastq_sbatch(sp_dir,config_info["sample_local_dict"])
    
    #Combine sbatch filenames into single object:
    sra_dl_sbatch_filenames = ncbi_sra_dl_sbatch_filenames + ena_sra_dl_sbatch_filenames + local_fastq_process_sbatch_filenames
     
    #Submit SRA read sbatch files, only allow 20 SRA jobs to run (or pend) at a time (set max_jobs)  
    max_jobs = 20
    sra_dl_jobids = []
    completed_jobids = {}
    job_count = 0
    #First submit up to the maximum number of jobs quickly
    if len(sra_dl_sbatch_filenames) > max_jobs:
        for i in range(0,max_jobs):
            sra_dl_jobids.append(sbatch_submit(sra_dl_sbatch_filenames[job_count]))
            job_count += 1
            sleep(1)
    else:
        for i in range(0,len(sra_dl_sbatch_filenames)):
            sra_dl_jobids.append(sbatch_submit(sra_dl_sbatch_filenames[job_count]))
            job_count += 1
            sleep(1)
    #Add an extra sleep to give sacct a chance to catch up
    sleep(20)
    #Then, enter while loop that will continue until the number of completed jobs matches the. number of sbatch files
    while len(completed_jobids) < len(sra_dl_sbatch_filenames):
        num_running = num_pend_run(sra_dl_jobids,start_date)
        while num_running < max_jobs and job_count < (len(sra_dl_sbatch_filenames)):
            sra_dl_jobids.append(sbatch_submit(sra_dl_sbatch_filenames[job_count]))
            job_count += 1
            sleep(20)
            num_running = num_pend_run(sra_dl_jobids,start_date)
        job_statuses = all_jobs_status(start_date)
        for job in sra_dl_jobids:
            if job not in completed_jobids:
                try:
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        print("Job %s completed"%job)
                except:
                    pass
        sleep(30)
    
    #After all jobs have finished, report which jobs failed
    for job in completed_jobids:
        if completed_jobids[job] != "COMPLETED":
            print("SRA download and parse job or local fastq process job %s failed with code: %s"%(job,completed_jobids[job]))
    
    #####Final processing and mapping   
    #Trim fastq files with NGmerge
    #Set Read Group information
    #Map fastq files to genome with BWA
    #Sort (convert to BAM), Index with GATK4 (installed in path)
    #Calculate alignment stats with GATK4 (installed in path)
    #Remove SAM file if stats file exists
    #Before submitting jobs, check to see that fastq files are there. If not, print info statement.
    mapping_jobids = []
    mapping_completed_jobids = {}
    
    for sample in config_info["sample_dict"]:
        for sra in config_info["sample_dict"][sample]:
            if os.path.isfile('%s/%s_1.fastq.gz'%(fastq_dir,sra)):
                if os.path.isfile('%s/%s_2.fastq.gz'%(fastq_dir,sra)):
                    if  os.path.isfile('%s/alignment/%s.sorted.bai'%(sp_dir,sra)) is False:
                        if os.path.isfile("%s/dedup/%s.dedup.sorted.bai"%(sp_dir,sample)) is False:
                            sra_map_sbatchfile = fastq_trim_align_stats(sp_dir,sra,config_info["abbv"],sample)
                            sra_map_jobid = sbatch_submit(sra_map_sbatchfile)
                            mapping_jobids.append(sra_map_jobid)
                            sleep(1)
                        else:
                            print('%s.dedup.sorted.bai already present, skipping'%(sample))
                    else:
                        print('%s.sorted.bai already present, skipping'%(sra))
                else:
                    print("No fastq_2 file for %s"%sra)
            else:
                print("No fastq files for %s"%sra)
    sleep(60)
       
    while len(mapping_completed_jobids) < len(mapping_jobids):
        job_statuses = all_jobs_status(start_date)
        for job in mapping_jobids:
            if job not in mapping_completed_jobids:
                if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                    mapping_completed_jobids[job] = job_statuses[job]
                    print("Job %s completed"%job)
        sleep(60)

    #After all jobs have finished, report which jobs failed
    for job in mapping_completed_jobids:
        if mapping_completed_jobids[job] != "COMPLETED":
            print("SRA trimming, mapping, sorting, and stats job %s failed with code: %s"%(job,mapping_completed_jobids[job]))
    

    
    now = datetime.datetime.now()
    print('Scripted finished: %s'%now)     
    
if __name__ == "__main__":
    main()
