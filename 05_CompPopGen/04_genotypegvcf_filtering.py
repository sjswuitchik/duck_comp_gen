#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
import shutil
from time import sleep
import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg') #Added to get working on the cluster
from matplotlib.backends.backend_pdf import PdfPages

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename, heterozygosity_arg, combine_gvcf_program_arg, bypass_interval_arg, memory_gg_arg, time_gg_arg):
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
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
            elif line[0] == "--HETEROZYGOSITY":
                if heterozygosity_arg is None:
                    config_info["het"] = line[1]
                else:
                    config_info["het"] = heterozygosity_arg
            elif line[0] == "--PIPELINE":
                config_info["pipeline"] = line[1]
            elif line[0] == "--NINTERVALS":
                config_info["nintervals"] = line[1]
            elif line[0] == "--MEMORY_HC":
                config_info["memory_hc"] = line[1]
            elif line[0] == "--TIME_HC":
                config_info["time_hc"] = line[1]
            elif line[0] == "--MEMORY_GG":
                if memory_gg_arg is None:
                    config_info["memory_gg"] = line[1]
                else:
                    config_info["memory_gg"] = memory_gg_arg
            elif line[0] == "--TIME_GG":
                if time_gg_arg is None:
                    config_info["time_gg"] = line[1]
                else:
                    config_info["time_gg"] = time_gg_arg
            elif line[0] == "--COMBINE_GVCF_PROGRAM":
                if combine_gvcf_program_arg is None:
                    config_info["combine_gvcf_program"] = line[1]
                else:
                    config_info["combine_gvcf_program"] = combine_gvcf_program_arg
            elif line[0] == "--BYPASS_INTERVAL":
                if bypass_interval_arg is None:
                    config_info["bypass_interval"] = line[1]
                else:
                    config_info["bypass_interval"] = bypass_interval_arg
                
    config_file.close()
    
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
        
    if "het" not in config_info:
        if heterozygosity_arg is None:
            config_info["het"] = "0.001"
            print("No heterozygosity specified, using default human value (0.001)")
        else:
            config_info["het"] = heterozygosity_arg
    
    if "pipeline" not in config_info:
        config_info["pipeline"] = "lowcoverage"
        #print("No pipeline specified (highcoverage or lowcoverage), using lowcoverage.")
        
    if config_info["pipeline"] != "highcoverage" and config_info["pipeline"] != "lowcoverage":
        sys.exit("Pipeline must be set to either 'highcoverage' or 'lowcoverage")
    
    if "nintervals" not in config_info:
        config_info["nintervals"] = "10"
        print("No specification for the number of intervals to analyze, using 10")
    
    if "memory_hc" not in config_info:
        config_info["memory_hc"] = "8"
        #print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")
    
    if "time_hc" not in config_info:
        config_info["time_hc"] = "12"
        #print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")
    
    if "memory_gg" not in config_info:
        if memory_gg_arg is None:
            config_info["memory_gg"] = "16"
            print("No specification of how much memory to use for GenotypeGVCF, using 16 GB by default")
        else:
            config_info["memory_gg"] = memory_gg_arg
    
    if "time_gg" not in config_info:
        if time_gg_arg is None:
            config_info["time_gg"] = "24"
            print("No specification of how much time to use for GenotypeGVCF, using 24 hours by default")
        else:
            config_info["time_gg"]= time_gg_arg
        
    if "combine_gvcf_program" not in config_info:
        if combine_gvcf_program_arg is None:
            config_info["combine_gvcf_program"] = "GenomicsDBImport"
            print("No specification of which program to use to combine gvcf files, using GenomicsDBImport by default")
        else:
            config_info["combine_gvcf_program"] = combine_gvcf_program_arg
        
    if "bypass_interval" not in config_info:
        if bypass_interval_arg is None:
            config_info["bypass_interval"] = "FALSE"
            print("BYPASS_INTERVAL set to FALSE, will require gvcfs from all samples for all intervals to proceed")
        else:
            config_info["bypass_interval"] = bypass_interval_arg
    
    
    #Return objects    
    return(config_info)


#Check if a directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):

        os.makedirs(dir)

#Create generic slurm script
def script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%j.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)
    
#Create generic slurm script for arrays (-o and -e for arrays)
def array_script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%A_%a.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%A_%a.err\n\n{cmd}
    '''
    return(slurm_script)

#Submit filename to slurm with sbatch with a given amount of time and memroy, returns job id number
def sbatch_submit(filename,memory,timelimit):
    proc = Popen('sbatch --mem %s --time %s:00:00 %s '%(memory,timelimit,filename),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    stdout = stdout.decode("utf-8","ignore")
    stdout = stdout.strip()
    stdout = stdout.strip('Submitted batch job ')
    return(stdout)

#Submit filename to slurm with sbatch with a given amount of time and memroy, returns job id number - includes first argument to include memory limit in java opts (use memory - 2)
def sbatch_submit_array(filename,memory,timelimit, array_nums):
    proc = Popen('sbatch --mem %s000 --time %s:00:00 --array=%s %s %d'%(memory,timelimit,array_nums,filename,(int(memory)-2)),shell=True,stdout=PIPE,stderr=PIPE)
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
    
#Count the number of pending jobs
def num_pend_run(job_id_list,date):
    count = 0
    status_dict = all_jobs_status(date)
    for job in job_id_list:
        if status_dict[job] == "PENDING" or status_dict[job] == "RUNNING":
            count += 1
    return(count)

#Check for missing gvcf interval files for a given sample in list of files. Will return a list of the missing intervals
def check_missing_gvcfs(arraystart,arrayend,sample_files,sample):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s.%s.g.vcf.gz"%(sample,str(i)) not in sample_files or "%s.%s.g.vcf.gz.tbi"%(sample,str(i)) not in sample_files:
            missing_ints.append(str(i))            
    
    return(missing_ints)
    
#Check for missing vcf interval files for a given sample in list of files. Will return a list of the missing intervals
def check_missing_vcfs(arraystart,arrayend,vcf_files,sp_abbr):
    missing_ints = []
    #Check if file_ext is a string, if so, just test that one type
    for i in range(arraystart,arrayend+1):
        if "%s_hardfilters.%s.vcf.gz"%(sp_abbr,str(i)) not in vcf_files or "%s_hardfilters.%s.vcf.gz.tbi"%(sp_abbr,str(i)) not in vcf_files:
            missing_ints.append(str(i))            
    return(missing_ints)

#Count number of intervals in split genome file
def count_intervals(nintervals,outputdir):
    if nintervals != "CHROMOSOME":
        nintervalfiles = nintervals
    else:
        filelist = os.listdir(outputdir)
        nintervalfiles = len(filelist)
    return(nintervalfiles)



###Create sbatch scripts
#Create a genotypegvcf sbatch file for a sample
def genotypegvcf_sbatch(sp_dir,sp_abbr,sample_list,het,nintervals,memory_gg,combine_gvcf_program):
    slurm_script = array_script_create()
    nintervals = str(nintervals)
    
    #Need to create string of --variant sample1.g.vcf.gz --variant sample2.g.vcf.gz for each sample in the list
    new_sample_list = []
    for i in range(0,len(sample_list)):
        new_sample_list.append('--variant %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz'%(sp_dir,sample_list[i]))

    all_sample_variant_call = " ".join(new_sample_list)

    
    #Load modules and get versions for all programs used
    ##For now, using my own installation of GATK as it is not yet installed on the cluster
    cmd_1 = 'module load jdk/1.8.0_45-fasrc01'
    
    cmd_2 = 'MEM=$1'    
    
    #Before running GenotypeGVCFs, need to combine GVCFs for all individuals into a single file
    
    if combine_gvcf_program == "CombineGVCFs":
        cmd_3 = 'if [ -f %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz.tbi ]\nthen\necho "gvcf file already complete, will not recreate"\nelse\nif [ -f %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz ]\nthen\nrm %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz\nelse\necho "No Combined gvcf file yet"\nfi\n/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" CombineGVCFs -R %s/genome/%s.fa %s -O %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list\nfi'%(sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,sp_abbr,all_sample_variant_call,sp_dir,sp_abbr,sp_dir,nintervals,sp_abbr)
    
        cmd_4 = '/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenotypeGVCFs -R %s/genome/%s.fa -V %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -O %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz --heterozygosity %s --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    elif combine_gvcf_program == "GenomicsDBImport":
        cmd_3 = '/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenomicsDBImport -R %s/genome/%s.fa %s --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list --genomicsdb-workspace-path %s/genomics_db/interval_${SLURM_ARRAY_TASK_ID}'%(sp_dir,sp_abbr,all_sample_variant_call,sp_dir,nintervals,sp_abbr,sp_dir)
    
        cmd_4 = '/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" GenotypeGVCFs -R %s/genome/%s.fa -V gendb://%s/genomics_db/interval_${SLURM_ARRAY_TASK_ID} -O %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz --heterozygosity %s --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,sp_dir,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    
    else:
        print("Program to combine gvcf files not specified")
        sys.exit()
    
    #Extract stat distributions    
    cmd_5 = '/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" VariantsToTable -V %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz -O %s/stats/%s_${SLURM_ARRAY_TASK_ID}_unfilteredVCFstats.txt -F CHROM -F POS -F TYPE -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F NCALLED -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -R %s/genome/%s.fa --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,nintervals,sp_abbr)
    
    #Apply default GATK hard filters
    cmd_6 = r"""/n/holyscratch01/informatics/swuitchik/CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" VariantFiltration -R %s/genome/%s.fa -V  %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz -O  %s/vcf/%s_hardfilters.${SLURM_ARRAY_TASK_ID}.vcf.gz --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0) " --filter-name "GATK_default" --filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" --filter-name "GATK_default" --filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" --filter-name "GATK_default" """%(sp_dir,sp_abbr,sp_dir,sp_abbr,sp_dir,sp_abbr)
    
    #Calculate missingness per individual
    cmd_7 = 'vcftools --gzvcf %s/vcf/%s_hardfilters.${SLURM_ARRAY_TASK_ID}.vcf.gz --missing-indv --out %s/stats/%s_ind_missingness.${SLURM_ARRAY_TASK_ID}'%(sp_dir,sp_abbr,sp_dir,sp_abbr)
    
    cmd_8 = 'rm %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz\nrm %s/vcf/%s.${SLURM_ARRAY_TASK_ID}.vcf.gz.tbi'%(sp_dir,sp_abbr,sp_dir,sp_abbr)
        
    cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7,cmd_8]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    genotypegvcf_script = slurm_script.format(partition="shared",cores="2",nodes="1",jobid="gg",sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/06_genotypegvcf_array.sbatch"%(sp_dir)
    out_file = open(out_filename,"w")
    out_file.write(genotypegvcf_script)
    out_file.close

    return(out_filename)
    


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    parser.add_argument("--HETEROZYGOSITY", help="heterozygosity to use with GenotypeGVCF, default 0.001. Setting value here will overwrite any values set in config file.", required=False)
    parser.add_argument("--COMBINE_GVCF_PROGRAM", help="use GenomicsDBImport or CombineGVCFs to combine gvcfs prior to genotyping, default GenomicsDBImport. Setting value here will overwrite any values set in config file", required=False)
    parser.add_argument("--BYPASS_INTERVAL", help="should missing gvcfs be allowed to continue with GenotypeGVCF? If TRUE will continue despite missing data, if FALSE will print a statement and exit the program (Default FALSE). Setting value here will overwrite any values set in config file", required=False)
    parser.add_argument("--MEMORY_GG", help="memory in GB (e.g. 8 = 8GB) to use for GenotypeGVCFs sbatch script, default 16. Setting value here will overwrite any values set in config file", required=False)
    parser.add_argument("--TIME_GG", help="time in hours (e.g. 12 = 12 hours) to use for GenotypeGVCFs sbatch script, default 24. Setting value here will overwrite any values set in config file", required=False)
    
    args = parser.parse_args()
    config_filename = args.config
    
    if args.HETEROZYGOSITY:
        heterozygosity_arg = args.HETEROZYGOSITY
        print("Using command line specified heterozygosity of %s"%heterozygosity_arg)
    else:
        heterozygosity_arg = None
    if args.COMBINE_GVCF_PROGRAM:
        combine_gvcf_program_arg = args.COMBINE_GVCF_PROGRAM
        print("Using command line specified COMBINE_GVCF_PROGRAM of %s"%combine_gvcf_program_arg)
    else:
        combine_gvcf_program_arg = None
    if args.BYPASS_INTERVAL:
        bypass_interval_arg = args.BYPASS_INTERVAL
        print("Using command line specified BYPASS_INTERVAL argument of %s"%(bypass_interval_arg))
    else:
        bypass_interval_arg = None
    if args.MEMORY_GG:
        memory_gg_arg = args.MEMORY_GG
        print("Using command line specified memory of %s GB for GenotypeGVCFs sbatch script"%memory_gg_arg)
    else:
        memory_gg_arg = None
    if args.TIME_GG:
        time_gg_arg = args.TIME_GG
        print("Using command line specified time of %s hours for GenotypeGVCFs sbatch script"%time_gg_arg)
    else:
        time_gg_arg = None
            
    now = datetime.datetime.now()
    print('Staring work on script 04: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get attributes
    config_info = extract_config(config_filename, heterozygosity_arg, combine_gvcf_program_arg, bypass_interval_arg, memory_gg_arg, time_gg_arg)

    #####Check if working directories exist, if not creates them
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    print("\nOutput will be written to %s\n"%sp_dir)
    
    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    dedup_dir = "%s/dedup"%(sp_dir)
    gvcf_dir = "%s/gvcf"%(sp_dir)
    vcf_dir = "%s/vcf"%(sp_dir)
    genomics_db_dir = "%s/genomics_db"%(sp_dir)
       
    directory_create(sp_dir)    
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(genome_dir)
    directory_create(stats_dir)
    directory_create(dedup_dir)
    directory_create(gvcf_dir)
    directory_create(vcf_dir)
    directory_create(genomics_db_dir)


    #Need to count intervals from genome split in previous step. If the genome was split by a number of intervals, this is trivial (it is just the number), but if it was split by chromosome, need to count files.
        
    nintervalfiles = count_intervals(config_info["nintervals"],"%s/%s_splits_interval_lists/"%(genome_dir,config_info["nintervals"]))
    
    
    #First, iterate through and make sure all HaplotypeCaller result files are present
    missing_dict = {}
    all_gvcfs = os.listdir(gvcf_dir)
    for sample in config_info["sample_dict"]:
        missing_ints = check_missing_gvcfs(1,int(nintervalfiles),all_gvcfs,sample)
        if len(missing_ints) > 0:
            missing_dict[sample] = missing_ints
        else:
            pass
    
    #If there are any missing files, print statement and kill script
    if len(missing_dict) > 0:
        print('There are missing gvcf files for given samples and interval set:')
        for sample in missing_dict:
            missing = ", ".join(missing_dict[sample])
            print('%s: %s'%(sample,missing))
        if config_info["bypass_interval"] == 'TRUE':
            pass
        else:
            sys.exit("Rerun missing intervals, change input samples, or set --BYPASS_INTERVAL TRUE in config file")
    else:
        pass
  

    #####Run GenotypeGVCF
    
    #Create GenotypeGVCF file - here we include all samples, so only a single slurm array is needed
    gg_filename = genotypegvcf_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample_list=list(config_info["sample_dict"].keys()),het=config_info["het"],nintervals=config_info["nintervals"],memory_gg=config_info["memory_gg"],combine_gvcf_program=config_info["combine_gvcf_program"])
    
    #Get number of finished files
    vcf_files = os.listdir(vcf_dir)
    finished_files = len([name for name in vcf_files if ".tbi" in name])
    
    all_jobids = []
    #Submit file the first time
    if finished_files == 0:
        
        #Submit job
        base_jobid = sbatch_submit_array(gg_filename,memory=config_info["memory_gg"],timelimit=config_info["time_gg"],array_nums="1-%s"%str(nintervalfiles))
        
        #Expand job IDs
        for i in range(1,int(nintervalfiles)+1):
            all_jobids.append("%s_%d"%(base_jobid,i))
    
    elif finished_files < int(nintervalfiles):
        #Check each interval, see if it has both a .vcf.gz and .tbi file
        
        missing = check_missing_vcfs(arraystart=1,arrayend=int(nintervalfiles),vcf_files=vcf_files,sp_abbr=config_info["abbv"])
        missing_vec = ",".join(missing)
        
         #Submit job, get base jobid for array
        base_jobid = sbatch_submit_array(gg_filename,memory=config_info["memory_gg"],timelimit=config_info["time_gg"], array_nums=missing_vec)
        sleep(1)
    
        #Add jobids for array to dictionary with jobid as key and sample as value
        for i in missing:
            all_jobids.append("%s_%s"%(base_jobid,i))
            
    elif finished_files == int(nintervalfiles):
        print("All vcf files present")
        
    else:
        sys.exit("More vcf files present than expected, check")
        
    #Give sacct a chance to catch up       
    sleep(20)
    
    #Then, enter while loop that will continue until the number of completed jobs matches the number of intervals
    #Create dictionary of completed jobids and completion statuses
    
    if finished_files != int(nintervalfiles):
        completed_jobids = {}
        rerun_jobids = []
        successful_intervals = []
        failed_intervals = []
    
        while len(completed_jobids) < len(all_jobids):
            job_statuses = all_jobs_status(start_date)
            current_jobs = all_jobids
            for job in current_jobs:
                if job not in completed_jobids:
                    if job in job_statuses:#Have to add this because array jobs may be delayed
                        if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                            completed_jobids[job] = job_statuses[job]
                            array_id = job.split("_")[1]
                        
                            #If job_id is "COMPLETED", check to make sure both the .vcf.gz file and .tbi file are both present. If they are, print and add to successful_samples dictionary (sample:[intervals])
                            if job_statuses[job] == "COMPLETED":
                                if os.path.isfile("%s/vcf/%s.%s.vcf.gz"%(sp_dir,config_info["abbv"],array_id)) and os.path.isfile("%s/vcf/%s.%s.vcf.gz.tbi"%(sp_dir,config_info["abbv"],array_id)):
                                    print("GenotypeGVCF job %s completed for interval %s"%(job, array_id))
                                    successful_intervals.append(array_id)

                            #If job_id is not COMPLETED, it means there was some sort of failure in the job. Resubmit with 2x time (up to 7 days, or 168 hours) and 2x memory
                            elif job_statuses[job] != "COMPLETED" and job not in rerun_jobids:
                                new_mem = str(int(config_info["memory_gg"])*2)
                                new_time =  int(config_info["time_gg"])*2
                                if new_time > 168:
                                   new_time = '168'
                                else:
                                    new_time = str(new_time)
                            
                                #Remove previous GenomicsDB if it exists - will throw error if directory already exists.
                                if os.path.isdir('%s/genomics_db/interval_%s'%(sp_dir,str(i))):
                                    shutil.rmtree(path = '%s/genomics_db/interval_%s'%(sp_dir,str(i)))
                            
                                #Submit array with only that interval
                                resubmitted_jobid = sbatch_submit_array(gg_filename,memory=new_mem,timelimit=new_time, array_nums=array_id)
                                sleep(1)
                            
                                #Add job id (including array number) to both rerun_jobids and all_jobids
                                rerun_jobids.append('%s_%s'%(resubmitted_jobid,array_id))
                            
                                all_jobids.append('%s_%s'%(resubmitted_jobid,array_id))
                            
                                print("GenotypeGVCF job %s failed, retrying interval %s with %s memory and %s time"%(job,array_id,new_mem,new_time))
                        
                            #If just doesn't finish and already resubmitted, do not submit again, print failure to log file, and add to failed_intervals list
                            elif job_statuses[job] != "COMPLETED" and job in rerun_jobids:
                                print("GenotypeGVCF job %s failure 2x for and interval %s"%(job,array_id))
                                failed_intervals.append(array_id)
                                
                            else:
                                print("Error with GenotypeGVCF job checking and resubmissions")
                        
            sleep(30)
    
        #After all jobs have finished, report intervals failed twice
        failed_intervals = ",".join(failed_intervals)
        if len(failed_intervals) > 0:
            print("Jobs failed for intervals: %s"%(failed_intervals))
    
    #Concatenate all missingness information into single file, adding additional column for easy manipulation in R, and additional file with mean and SD missingness per individual
    all_missing_file = open('%s/stats/_%s_all_all_missingness_info.txt'%(sp_dir,config_info["abbv"]),'w')
    mean_sd_missing_file = open('%s/stats/_%s_all_mean_missingness_info.txt'%(sp_dir,config_info["abbv"]),'w')
    
    sample_miss_dict = {}
    
    #Write headers
    all_missing_file.write('INDV	N_DATA	N_GENOTYPES_FILTERED	N_MISS	F_MISS	INTERVAL\n')
    mean_sd_missing_file.write('INDV	MEAN_MISS	SD_MISS	MAX_MISS	MIN_MISS\n')
    
    #If missingness file exists for interval, write results to all_missing_file, and add fraction missing to dictionary for each sample
    for i in range(1,int(nintervalfiles)+1):
        
        if os.path.isfile('%s/stats/%s_ind_missingness.%d.imiss'%(sp_dir,config_info["abbv"],i)):
            missing_file = open('%s/stats/%s_ind_missingness.%d.imiss'%(sp_dir,config_info["abbv"],i),'r')
            for line in missing_file:
                line = line.strip()
                split_line = line.split()
                if split_line[0] != 'INDV':
                    all_missing_file.write('%s\t%d\n'%(line,i))
                    if split_line[0] in sample_miss_dict:
                        sample_miss_dict[split_line[0]].append(float(split_line[4]))
                    else:
                        sample_miss_dict[split_line[0]] = [float(split_line[4])]
            
            missing_file.close()
    
    for sample in sample_miss_dict:
        sample_mean = np.round(np.mean(sample_miss_dict[sample]),3)
        sample_sd = np.round(np.std(sample_miss_dict[sample]),3)
        sample_min = np.min(sample_miss_dict[sample])
        sample_max = np.max(sample_miss_dict[sample])
        
        mean_sd_missing_file.write('%s\t%s\t%s\t%s\t%s\n'%(sample,str(sample_mean),str(sample_sd),str(sample_max),str(sample_min)))
    
    all_missing_file.close()
    mean_sd_missing_file.close()
    
    #Copy summary missingness file over to summary directory
    general_dir = "_ALL_SPECIES_SUMMARIES"
    
    try:
        proc = Popen('cp %s/stats/_%s_all_mean_missingness_info.txt %s/%s_all_mean_missingness_info.txt'%(sp_dir,config_info["abbv"],general_dir,config_info["abbv"]),shell=True)
    except:
        print("There was an error copying summary missingness file")
        
    #For each interval, check to make sure vcf .tbi file exists, and if so, delete either combined gvcf (if program is CombineGVCFs) or the GenomicsDB (if GenomicsDBImport was used). If no .tbi file exists, prints statement to log file.    
    for interval in range(1,int(nintervalfiles)+1):
        if os.path.isfile("%s/vcf/%s_hardfilters.%d.vcf.gz.tbi"%(sp_dir,config_info["abbv"],interval)):
            if config_info["combine_gvcf_program"] == "CombineGVCFs":
                if os.path.isfile("%s/gvcf/%s.%d.gvcf.gz.tbi"%(sp_dir,config_info["abbv"],interval)) and os.path.isfile("%s/gvcf/%s.%d.gvcf.gz"%(sp_dir,config_info["abbv"],interval)):
                    proc = Popen("rm %s/gvcf/%s.%d.gvcf.gz.tbi"%(sp_dir,config_info["abbv"],interval),shell=True)
                    proc = Popen("rm %s/gvcf/%s.%d.gvcf.gz"%(sp_dir,config_info["abbv"],interval),shell=True)
            elif config_info["combine_gvcf_program"] == "GenomicsDBImport":
                if os.path.isdir('%s/genomics_db/interval_%d'%(sp_dir,interval)):
                    shutil.rmtree(path = '%s/genomics_db/interval_%d'%(sp_dir,interval))
        else:
            print("No vcf .tbi file for interval %d, check"%(interval))
            
                
    now = datetime.datetime.now()
    print('Finished script 04: %s'%now)


if __name__ == "__main__":
    main()

