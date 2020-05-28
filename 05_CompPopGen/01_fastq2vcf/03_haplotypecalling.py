#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg') #Added to get working on the cluster
from matplotlib.backends.backend_pdf import PdfPages

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename, heterozygosity_arg, pipeline_arg, memory_hc_arg, time_hc_arg):
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
                if pipeline_arg is None:
                    config_info["pipeline"] = line[1]
                else:
                    config_info["pipeline"] = pipeline_arg
            elif line[0] == "--NINTERVALS":
                config_info["nintervals"] = line[1]
            elif line[0] == "--MEMORY_HC":
                if memory_hc_arg is None:
                    config_info["memory_hc"] = line[1]
                else:
                    config_info["memory_hc"] = memory_hc_arg
            elif line[0] == "--TIME_HC":
                if time_hc_arg is None:
                    config_info["time_hc"] = line[1]
                else:
                    config_info["time_hc"] = time_hc_arg
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
        if pipeline_arg is None:
            config_info["pipeline"] = "lowcoverage"
            print("No pipeline specified (highcoverage or lowcoverage), using lowcoverage.")
        else:
            config_info["pipeline"] = pipeline_arg
        
    if config_info["pipeline"] != "highcoverage" and config_info["pipeline"] != "lowcoverage":
        sys.exit("Pipeline must be set to either 'highcoverage' or 'lowcoverage'")
    
    if "nintervals" not in config_info:
        config_info["nintervals"] = "10"
        print("No specification for the number of intervals to analyze, using 10")
    
    if "memory_hc" not in config_info:
        if memory_hc_arg is None:
            config_info["memory_hc"] = "8"
            print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")
        else:
            config_info["memory_hc"] = memory_hc_arg
    
    if "time_hc" not in config_info:
        if time_hc_arg is None:
            config_info["time_hc"] = "12"
            print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")
        else:
            config_info["time_hc"] = time_hc_arg
    
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

##Currently requires NCBI names if you specify chromosome. Add config option to specify regex to identify chromosomes, have default be NCBI.
def split_genome(sp_dir,sp_abbr,nintervals,outputdir):
	#Open input		
    fai = open("%s/genome/%s.fa.fai"%(sp_dir,sp_abbr),"r")
    outHandle = sp_abbr

    if nintervals != "CHROMOSOME":
        faiList = []
        totLen = 1
        cumStart = []
        cumEnd = []
    
        n = int(nintervals)
    
        #For each line in the .fai file, append data to a list, and create two additional lists with cumulative start and end positions.
        for line in fai:
            line = line.strip().split("\t")
            faiList.append(line)
            cumStart.append(totLen)
            totLen = totLen + int(line[1])
            cumEnd.append(totLen)

        #Create interval based on total length of scaffolds, and create list of start and end values
        #### Would be better to refactor to just compute interval size, and add scaffolds to intervals until the interval exceeds size, then move on to the next interval.
        interval = (totLen/n)+1
        intervalNums = range(1,n+1)
        intervalEnd = [x * interval for x in intervalNums]
        intervalStart = [(x - interval) + 1 for x in intervalEnd]
    
        #Create list of which file each scaffold should go to, depending on where it falls on the interval list. Note that scaffolds that span intervals will be moved to the preceeding interval.
        fileDes = []
        for i in range(len(faiList)):
            fileDes.append(0)
            for j in range(n):
                if cumStart[i] >= intervalStart[j] and cumEnd[i] < intervalEnd[j]:
                    fileDes[i] = (j+1)
                elif j < (n-1) and cumStart[i] >= intervalStart[j] and cumEnd[i] < intervalEnd[(j+1)] and cumStart[i] < intervalStart[(j+1)]:
                    fileDes[i] = (j+1)
    
        #Create list of output files based on species abbreviation
        outList = []
    
        for i in range(n):
            outFile = (outputdir+outHandle+"_"+(str(i+1))+".interval_list")
            outList.append(open(outFile,"w"))
        
        #Append scaffold name to appropriate output file.						
        for i in range(len(faiList)):
            outList[fileDes[(i)]-1].write(faiList[i][0]+"\n")
        
        for i in range(n):
            outList[i].close()
        
        nintervalfiles = int(nintervals)
    
    else:
        faiList = []
        for line in fai:
            line = line.strip().split("\t")
            faiList.append(line)
        #Get number of chromosomes, will produce that many interval files + 1
        nchr = len([name for name in faiList if ("NC_" in name[0] or "CM" in name[0])])
        nintervalfiles=(nchr+1)
        filenum = 1
        #Open file to write all non-chromosomes
        randoutFile = open((outputdir+outHandle+"_"+(str(nintervalfiles))+".interval_list"),"w")
        for i in range(len(faiList)):
            if "NC_" in faiList[i][0] or "CM" in faiList[i][0]:
                outFile = open((outputdir+outHandle+"_"+(str(filenum))+".interval_list"),"w")
                outFile.write(faiList[i][0])
                outFile.close()
                filenum += 1
            else:
                randoutFile.write("%s\n"%faiList[i][0])
                
        randoutFile.close()

    fai.close()
	
    return(nintervalfiles)



###Create sbatch scripts
#Create a haplotypecaller sbatch file for a sample
def haplotypecaller_sbatch(sp_dir,sp_abbr,sample,het,memory_hc,nintervals,pipeline):
    slurm_script = array_script_create()
    nintervals = str(nintervals)

    #Load modules and get versions for all programs used
    ##For now, using my own installation of GATK as it is not yet installed on the cluster
    cmd_1 = 'module load jdk/1.8.0_45-fasrc01'
    
    cmd_2 = 'MEM=$1'
    
    if pipeline == "highcoverage":
    #Command to donwsample if proportion <0.95, if >0.95, just copy
        cmd_3 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" HaplotypeCaller -I %s/dedup/%s.dedup.sorted.bam -O %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -R %s/genome/%s.fa --heterozygosity %s --ERC GVCF --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list'%(sp_dir,sample,sp_dir,sample,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    
    elif pipeline == "lowcoverage":
        cmd_3 = 'gatk --java-options "-Xmx${MEM}g -XX:ParallelGCThreads=1" HaplotypeCaller -I %s/dedup/%s.dedup.sorted.bam -O %s/gvcf/%s.${SLURM_ARRAY_TASK_ID}.g.vcf.gz -R %s/genome/%s.fa --heterozygosity %s --ERC GVCF --intervals %s/genome/%s_splits_interval_lists/%s_${SLURM_ARRAY_TASK_ID}.interval_list --min-dangling-branch-length 1 --min-pruning 1'%(sp_dir,sample,sp_dir,sample,sp_dir,sp_abbr,het,sp_dir,nintervals,sp_abbr)
    
    cmd_list = [cmd_1,cmd_2,cmd_3]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    haplocaller_script = slurm_script.format(partition="shared",cores="2",nodes="1",jobid="hc",sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/05_haplotypecaller_%s_array.sbatch"%(sp_dir,sample)
    out_file = open(out_filename,"w")
    out_file.write(haplocaller_script)
    out_file.close

    return(out_filename)


    


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    parser.add_argument("--HETEROZYGOSITY", help="heterozygosity to use with HaplotypeCaller, default 0.001. Setting value here will overwrite any values set in config file.", required=False)
    parser.add_argument("--PIPELINE", help="highcoverage or lowcoverage set of HaplotypeCaller arguments, default lowcoverage. Setting value here will overwrite any values set in config file", required=False)
    parser.add_argument("--MEMORY_HC", help="Memory in GB (e.g. 8 = 8GB) to use for HaplotypeCaller sbatch script, default 16. Setting value here will overwrite any values set in config file", required=False)
    parser.add_argument("--TIME_HC", help="Time in hours (e.g. 12 = 12 hours) to use for HaplotypeCaller sbatch script, default 12. Setting value here will overwrite any values set in config file", required=False)
    args = parser.parse_args()
    config_filename = args.config
    if args.HETEROZYGOSITY:
        heterozygosity_arg = args.HETEROZYGOSITY
        print("Using command line specified heterozygosity of %s"%heterozygosity_arg)
    else:
        heterozygosity_arg = None
    if args.PIPELINE:
        pipeline_arg = args.PIPELINE
        print("Using command line specified pipeline of %s"%pipeline_arg)
    else:
        pipeline_arg = None
    if args.MEMORY_HC:
        memory_hc_arg = args.MEMORY_HC
        print("Using command line specified memory of %s GB for HaplotypeCaller sbatch script"%memory_hc_arg)
    else:
        memory_hc_arg = None
    if args.TIME_HC:
        time_hc_arg = args.TIME_HC
        print("Using command line specified time of %s hours for HaplotypeCaller sbatch script"%time_hc_arg)
    else:
        time_hc_arg = None
    
    now = datetime.datetime.now()
    print('Staring work on script 03: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get attributes
    config_info = extract_config(config_filename, heterozygosity_arg, pipeline_arg, memory_hc_arg, time_hc_arg)

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
       
    directory_create(sp_dir)    
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(genome_dir)
    directory_create(stats_dir)
    directory_create(dedup_dir)
    directory_create(gvcf_dir)
    directory_create(vcf_dir)


    #####Split genome into N specified intervals (if CHROMOSOME instead of number, will split into chromosomes, with all unplaced scaffolds in one file), create directory to store interval files. Function returns the number of interval files (nintervalfiles)
    #Added sleep to give proc time to copy the .fai file

    directory_create('%s/%s_splits_interval_lists/'%(genome_dir,config_info["nintervals"]))
    nintervalfiles = split_genome(sp_dir,config_info["abbv"],config_info["nintervals"],"%s/%s_splits_interval_lists/"%(genome_dir,config_info["nintervals"]))
  

    #####Run HaplotypeCaller
    
    #Submit all jobs the first time
    #hc_filenames is a dictionary with the sample as key and filename as value
    hc_filenames = {}
    #all_jobids is a dictionary with jobid (including array numbers as key and sample as value)
    all_jobids = {}

    for sample in config_info["sample_dict"]:

        sample_files = [name for name in os.listdir(gvcf_dir) if sample in name]
        finished_files = len([name for name in sample_files if ".tbi" in name])
        
        #If no files exist, submit full array
        if finished_files == 0:

            #Create sbatch file, add to filename dictionary with sample as key and filename as value
            if os.path.isfile('%s/dedup/%s.dedup.sorted.bai'%(sp_dir,sample)) and os.path.isfile('%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)):
                hc_filename = haplotypecaller_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample=sample,het=config_info["het"],memory_hc=config_info["memory_hc"],nintervals=config_info["nintervals"],pipeline=config_info["pipeline"])
                hc_filenames[sample] = hc_filename
            else:
                print('Deduped bam or bai file missing for sample %s, will not run HaplotypeCaller for that sample.'%(sample))
            
            #Submit job, get base jobid for array
            base_jobid = sbatch_submit_array(hc_filename,memory=config_info["memory_hc"],timelimit=config_info["time_hc"], array_nums="1-%d"%nintervalfiles)
            sleep(1)
        
            #Add jobids for array to dictionary with jobid as key and sample as value
            for i in range(1,nintervalfiles+1):
                all_jobids["%s_%d"%(base_jobid,i)] = sample
        
        #If the number of sample files is less than the the number of interval files x2 (because of vcf and index), that means some intervals are missing. Only submit those intervals that don't have .tbi (index) files.
        elif finished_files < int(nintervalfiles):
            #Check each interval, see if it has both a .vcf.gz and .tbi file
            hc_filename = haplotypecaller_sbatch(sp_dir,sp_abbr=config_info["abbv"],sample=sample,het=config_info["het"],memory_hc=config_info["memory_hc"],nintervals=config_info["nintervals"],pipeline=config_info["pipeline"])
            hc_filenames[sample] = hc_filename
            
            missing = check_missing_gvcfs(arraystart=1,arrayend=nintervalfiles,sample_files=sample_files,sample=sample)
            
            missing_vec = ",".join(missing)

        
             #Submit job, get base jobid for array
            base_jobid = sbatch_submit_array(hc_filename,memory=config_info["memory_hc"],timelimit=config_info["time_hc"], array_nums=missing_vec)
            sleep(1)
        
            #Add jobids for array to dictionary with jobid as key and sample as value
            for i in missing:
                all_jobids["%s_%s"%(base_jobid,i)] = sample
            
        elif finished_files == int(nintervalfiles):
            print("Sample %s has all gvcf files, skipping HaplotypeCaller"%sample)
        
        else:
            print("Sample %s has more gvcfs than expected, check"%sample)
            
    
    #Give sacct a chance to catch up       
    sleep(20)
    
    #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
    #Create dictionary of completed jobids and completion statuses
    completed_jobids = {}
    rerun_jobids = {}
    successful_samples = {}
    failed_samples = {}
    
    while len(completed_jobids) < len(all_jobids):
        job_statuses = all_jobs_status(start_date)
        current_jobs = list(all_jobids.keys())
        for job in current_jobs:
            if job not in completed_jobids:
                if job in job_statuses:#Have to add this because array jobs may be delayed
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        array_id = job.split("_")[1]
                        
                        #If job_id is "COMPLETED", check to make sure both the .vcf.gz file and .tbi file are both present. If they are, print and add to successful_samples dictionary (sample:[intervals])
                        if job_statuses[job] == "COMPLETED":
                            if os.path.isfile("%s/gvcf/%s.%s.g.vcf.gz"%(sp_dir,all_jobids[job],array_id)) and os.path.isfile("%s/gvcf/%s.%s.g.vcf.gz.tbi"%(sp_dir,all_jobids[job],array_id)):
                                print("Job %s completed for sample %s"%(job, all_jobids[job]))
                                if all_jobids[job] in successful_samples:
                                    successful_samples[all_jobids[job]].append(array_id)
                                else:
                                    successful_samples[all_jobids[job]] = [array_id]
                        #If job_id is not COMPLETED, it means there was some sort of failure in the job. Resubmit with 2x time (up to 7 days, or 168 hours) and 2x memory
                        elif job_statuses[job] != "COMPLETED" and job not in rerun_jobids:
                            new_mem = str(int(config_info["memory_hc"])*2)
                            new_time =  int(config_info["time_hc"])*2
                            if new_time > 168:
                               new_time = '168'
                            else:
                                new_time = str(new_time)
                            #Submit array with only that interval
                            resubmitted_jobid = sbatch_submit_array(hc_filenames[all_jobids[job]],memory=new_mem,timelimit=new_time, array_nums=array_id)
                            sleep(1)
                            
                            #Add job id (including array number) to both rerun_jobids and all_jobids
                            rerun_jobids['%s_%s'%(resubmitted_jobid,array_id)] = all_jobids[job]
                            
                            all_jobids['%s_%s'%(resubmitted_jobid,array_id)] = all_jobids[job]
                            
                            print("Job %s failed, retrying sample %s, interval %s with %s memory and %s time"%(job, all_jobids[job],array_id,new_mem,new_time))
                        
                        #If just doesn't finished and already resubmitted, do not submit again, print failure to log file, and add to failed_samples dictionary
                        elif job_statuses[job] != "COMPLETED" and job in rerun_jobids:
                            print("HaplotypeCaller job %s failure 2x for sample %s and interval %s"%(job,all_jobids[job],array_id))    
                            if all_jobids[job] in failed_samples:
                                failed_samples[all_jobids[job]].append(array_id)
                            else:
                                failed_samples[all_jobids[job]] = [array_id]
                                
                        else:
                            print("Error with HaplotypeCaller job checking and resubmissions")
                        
        sleep(30)
    
    #After all jobs have finished, report which samples and intervals failed twice
    for sample in failed_samples:
        failed_intervals = ",".join(failed_samples[sample])
        print("Sample %s, failed for intervals: %s"%(sample,failed_intervals))

               
    now = datetime.datetime.now()
    print('Finished script 03: %s'%now)


if __name__ == "__main__":
    main()

