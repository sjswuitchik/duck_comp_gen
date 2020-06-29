#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
import gzip
from subprocess import Popen,PIPE
from time import sleep
import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg') #Added to get working on the cluster
from matplotlib.backends.backend_pdf import PdfPages

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
            elif line[0] == "--OUT_DIR":
                config_info["out_dir"] = line[1]
            elif line[0] == "--HETEROZYGOSITY":
                config_info["het"] = line[1]
            elif line[0] == "--PIPELINE":
                config_info["pipeline"] = line[1]
            elif line[0] == "--NINTERVALS":
                config_info["nintervals"] = line[1]
            elif line[0] == "--MEMORY_HC":
                config_info["memory_hc"] = line[1]
            elif line[0] == "--TIME_HC":
                config_info["time_hc"] = line[1]
            elif line[0] == "--MEMORY_GG":
                config_info["memory_gg"] = line[1]
            elif line[0] == "--TIME_GG":
                config_info["time_gg"] = line[1]
            elif line[0] == "--COMBINE_GVCF_PROGRAM":
                config_info["combine_gvcf_program"] = line[1]
            elif line[0] == "--BYPASS_INTERVAL":
                config_info["bypass_interval"] = line[1]

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
        config_info["het"] = "0.001"
        #print("No heterozygosity specified, using default human value (0.001)")

    if "pipeline" not in config_info:
        config_info["pipeline"] = "lowcoverage"
        #print("No pipeline specified (highcoverage or lowcoverage), using lowcoverage.")

    #if config_info["pipeline"] != "highcoverage" and config_info["pipeline"] != "lowcoverage":
        #sys.exit("Pipeline must be set to either 'highcoverage' or 'lowcoverage")

    if "nintervals" not in config_info:
        config_info["nintervals"] = "10"
        #print("No specification for the number of intervals to analyze, using 10")

    if "memory_hc" not in config_info:
        config_info["memory_hc"] = "8"
        #print("No specification of how much memory to use for HaplotypeCaller, using 8GB by default")

    if "time_hc" not in config_info:
        config_info["time_hc"] = "12"
        #print("No specification of how much time to use for HaplotypeCaller, using 12 hours by default")

    if "memory_gg" not in config_info:
        config_info["memory_gg"] = "8"
        #print("No specification of how much memory to use for GenotypeGVCF, using 8 GB by default")

    if "time_gg" not in config_info:
        config_info["time_gg"] = "12"
        #print("No specification of how much time to use for GenotypeGVCF, using 12 hours by default")

    if "combine_gvcf_program" not in config_info:
        config_info["combine_gvcf_program"] = "GenomicsDBImport"
        #print("No specification of which program to use to combine gvcf files, using GenomicsDBImport by default")

    if "bypass_interval" not in config_info:
        config_info["bypass_interval"] = "FALSE"
        #print("BYPASS_INTERVAL set to FALSE, will require gvcfs from all samples for all intervals to proceed")


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
    proc = Popen('sbatch --mem %s000 --time %s:00:00 %s '%(memory,timelimit,filename),shell=True,stdout=PIPE,stderr=PIPE)
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
        if "%s.%s.vcf.gz"%(sp_abbr,str(i)) not in vcf_files or "%s.%s.vcf.gz.tbi"%(sp_abbr,str(i)) not in vcf_files:
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



#Create an sbatch file for a given set of samples to create genome coverage graphs with bedtools - will use deduplicated BAM file and will write results to stats_coverage

def sample_coverage_sbatch(sp_dir,sp_abbr,sample):
    slurm_script = script_create()

    #First check if dedup file is present, if it is, continue and if not print statement.
    dedup_sorted_filename = '%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)
    genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
    if os.path.isfile(dedup_sorted_filename):

        #Load modules and get versions for all programs used
        cmd_1 = 'module load bedtools2/2.26.0-fasrc01'

        #Create bedgraph with bedtools of coverage (include regions with 0 coverage)
        cmd_2 = 'bedtools genomecov -bga -ibam %s/dedup/%s.dedup.sorted.bam -g %s/genome/%s.fa > %s/stats_coverage/%s.bg'%(sp_dir,sample,sp_dir,sp_abbr,sp_dir,sample)


        cmd_list = [cmd_1,cmd_2]

        final_cmd = "\n\n".join(cmd_list)


#Format sbatch script
        sample_coverage_script = slurm_script.format(partition="shared",time="0-12:00",mem="8000",cores="1",nodes="1",jobid="genomecov",sp_dir=sp_dir,cmd=final_cmd)
        out_filename = "%s/scripts/07_coverage_bedgraph_%s.sbatch"%(sp_dir,sample)
        out_file = open(out_filename,"w")
        out_file.write(sample_coverage_script)
        out_file.close
        return(out_filename)

    else:
        print('No sorted dedup file for %s, cannot compute coverage bedgraph.'%(sample))



#Create an sbatch file for a given set of samples to create genome coverage graphs with bedtools - will use deduplicated BAM file and will write results to stats_coverage

def union_coverage_sbatch(sp_dir,sp_abbr,sample_bedgraph_file_list,sample_names_bedgraph_file_list):
    slurm_script = script_create()

    sample_bedgraph_files = " ".join(sample_bedgraph_file_list)
    sample_names_bedgraph_files = " ".join(sample_names_bedgraph_file_list)
    #First check if dedup file is present, if it is, continue and if not print statement.
            #Load modules and get versions for all programs used
    cmd_1 = 'module load bedtools2/2.26.0-fasrc01'

    #Create bedgraph with bedtools of coverage (include regions with 0 coverage)
    cmd_2 = 'bedtools unionbedg -header -empty -g %s/genome/%s.fa.fai -i %s -names %s > %s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir,sp_abbr,sample_bedgraph_files,sample_names_bedgraph_files,sp_dir,sp_abbr)

    cmd_3 = 'gzip %s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir,sp_abbr)

    cmd_list = [cmd_1,cmd_2,cmd_3]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    sample_coverage_script = slurm_script.format(partition="shared",cores="1",nodes="1",jobid="genomecov",sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/08_union_coverage_bedgraph_all_samples.sbatch"%(sp_dir)
    out_file = open(out_filename,"w")
    out_file.write(sample_coverage_script)
    out_file.close
    return(out_filename)


#Takes in species directory and abbreviation, calls supplemental python script that will be run on cluster  ######Note that this is using my own python environment, might want to make universal for others in future.
def sum_coverage_sbatch(sp_dir,sp_abbr):
    slurm_script = script_create()

    #Load modules and get versions for all programs used
    cmd_1 = 'module load bedtools2/2.26.0-fasrc01\nmodule load Anaconda/5.0.1-fasrc01\nsource activate pipeline'

    #run python script to sum, create histogram if missing, and create clean sites, too high sites, and too low sites bedfiles
    cmd_2 = 'python 05_01_sum_coverage_subscript.py --sp_dir %s --sp_abbr %s'%(sp_dir,sp_abbr)

    #sort and merge all bedgraphs
    cmd_3 = 'sort -k 1,1 -k2,2n %s/stats_coverage/_%s_clean_coverage_sites.bed | bedtools merge -i stdin > %s/stats_coverage/_%s_clean_coverage_sites_merged.bed'%(sp_dir,sp_abbr,sp_dir,sp_abbr)

    cmd_4 = 'sort -k 1,1 -k2,2n %s/stats_coverage/_%s_too_low_coverage_sites.bed | bedtools merge -i stdin > %s/stats_coverage/_%s_too_low_coverage_sites_merged.bed'%(sp_dir,sp_abbr,sp_dir,sp_abbr)

    cmd_5 = 'sort -k 1,1 -k2,2n %s/stats_coverage/_%s_too_high_coverage_sites.bed | bedtools merge -i stdin > %s/stats_coverage/_%s_too_high_coverage_sites_merged.bed'%(sp_dir,sp_abbr,sp_dir,sp_abbr)

    cmd_list = [cmd_1,cmd_2,cmd_3,cmd_4,cmd_5]

    final_cmd = "\n\n".join(cmd_list)

#Format sbatch script
    sample_coverage_script = slurm_script.format(partition="shared",cores="1",nodes="1",jobid="sumcov",sp_dir=sp_dir,cmd=final_cmd)
    out_filename = "%s/scripts/08_sum_coverage_clean_bedgraph.sbatch"%(sp_dir)
    out_file = open(out_filename,"w")
    out_file.write(sample_coverage_script)
    out_file.close
    return(out_filename)


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping",  required=True)
    args = parser.parse_args()
    config_filename = args.config

    now = datetime.datetime.now()
    print('Staring work on script 05: %s'%now)
    start_date = now.strftime("%Y-%m-%d")

    #Open config file and get attributes
    config_info = extract_config(config_filename)

    #####Check if working directories exist, if not creates them
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])

    print("\nOutput will be written to %s\n"%sp_dir)

    logs_dir = "%s/logs"%(sp_dir)
    scripts_dir = "%s/scripts"%(sp_dir)
    genome_dir = "%s/genome"%(sp_dir)
    stats_dir = "%s/stats"%(sp_dir)
    dedup_dir = "%s/dedup"%(sp_dir)
    stats_coverage_dir = "%s/stats_coverage"%(sp_dir)


    directory_create(sp_dir)
    directory_create(logs_dir)
    directory_create(scripts_dir)
    directory_create(stats_coverage_dir)



    #Submit all jobs
    #Coverage_filenames is a dictionary with the sample as key and filename as value
    coverage_filenames = {}
    #all_jobids is a dictionary with jobid (including array numbers as key and sample as value)
    all_jobids = {}
    failed_samples = []

    for sample in config_info["sample_dict"]:
        genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
        if os.path.isfile(genome_cov_filename) and os.path.getsize(genome_cov_filename) > 0:
            print("Coverage file for sample %s already exists, will not recalculate"%(sample))
        else:
            coverage_filenames[sample] = sample_coverage_sbatch(sp_dir,config_info["abbv"],sample)
    if len(coverage_filenames) > 0:
        #Submit dedup read sbatch files
        coverage_jobids = {}
        completed_jobids = {}
        for sample in coverage_filenames:
            all_jobids[sbatch_submit(coverage_filenames[sample],memory=8,timelimit=8)] = sample
            sleep(1)
        #Add an extra sleep to give sacct a chance to catch up
        sleep(20)
        #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
        while len(completed_jobids) < len(coverage_filenames):
            num_running = num_pend_run(all_jobids,start_date)
            job_statuses = all_jobs_status(start_date)
            for job in all_jobids:
                if job not in completed_jobids:
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        print("Job %s completed for sample %s"%(job,all_jobids[job]))
            sleep(30)

        #After all jobs have finished, report which jobs failed, also add samples that failed to list to exclude from whole-genome calculations
        for job in completed_jobids:
            if completed_jobids[job] != "COMPLETED":
                print("Coverage bedgraph job %s failed with code: %s for sample %s"%(job,completed_jobids[job],all_jobids[job]))
                failed_samples.append(all_jobids[job])


    print("Now computing whole-genome coverage bedgraph\n")

    ########## Create union bedgraph with all sample coverages included
    #Make a list of all sample bedgraphs to include in genome coverage calculations. Iterates through samples names, makes sure coverage bedgraph files exist (if not, skips sample), and that sample completed above command successfully.
    sample_bedgraph_file_list = []
    sample_names_bedgraph_file_list = []
    for sample in config_info["sample_dict"]:
        genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
        if sample in failed_samples:
            print("Will not include sample %s in genome coverage calculation, job failed")
        elif os.path.isfile(genome_cov_filename) and os.path.getsize(genome_cov_filename) > 0:
            sample_bedgraph_file_list.append(genome_cov_filename)
            sample_names_bedgraph_file_list.append(sample)

    #Checks if union bedgraph file exists - if so, skips creation.
    if os.path.isfile('%s/stats_coverage/_%s_all_samples_union.bg.gz'%(sp_dir,config_info["abbv"])) and os.path.getsize('%s/stats_coverage/_%s_all_samples_union.bg.gz'%(sp_dir,config_info["abbv"])) > 0:
        print("Union bedfile already exists, will not recreate")

    else:
        #if there is more than one sample, create and submit bedgraph job
        if len(sample_bedgraph_file_list) > 1:
            #Create and submit file for union bedgraph job
            union_sbatch_file = union_coverage_sbatch(sp_dir,config_info["abbv"],sample_bedgraph_file_list,sample_names_bedgraph_file_list)
            union_job_id = sbatch_submit(union_sbatch_file,memory=8,timelimit=72)
            sleep(30)

            #Only check on union job if actually submitted.
            if union_job_id is not None:
                dones = ['COMPLETED','CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL']
                #Check job id status of union job. If not in one of the 'done' job status categories, wait 30 seconds and check again.
                while jobid_status(union_job_id,start_date) not in dones:
                    sleep(30)

                #Check to make sure job completed, and that all necessary files are present. If not, exit and give information.
                union_job_completion_status = jobid_status(union_job_id,start_date)
                if union_job_completion_status != 'COMPLETED':
                    sys.exit("There was a problem creating the union bedgraph file. The job exited with status %s. Please diagnose and fix before moving on"%union_job_completion_status)

                #If job successfully completed, union gzipped file exists and has a size > 0, remove all sample bedgraphs if they are still present
                else:
                    if os.path.isfile('%s/stats_coverage/_%s_all_samples_union.bg.gz'%(sp_dir,config_info["abbv"])) and os.path.getsize('%s/stats_coverage/_%s_all_samples_union.bg.gz'%(sp_dir,config_info["abbv"])) > 0:
                        for sample in config_info["sample_dict"]:
                            genome_cov_filename = '%s/stats_coverage/%s.bg'%(sp_dir,sample)
                            if os.path.isfile(genome_cov_filename):
                                proc = Popen('rm %s'%genome_cov_filename,shell=True)
        else:
            #If only one sample, create file with header, then just cat the bedgraph from that one sample.
            #Create union header file for one sample
            union_header_file = open('%s/stats_coverage/_%s_all_samples_union_header.bg'%(sp_dir,config_info["abbv"]), "w")
            union_header_file.write('chrom\tstart\tend\t%s\n'%(sample_names_bedgraph_file_list[0]))
            union_header_file.close()

            #Cat header file with bedgraph of that sample
            proc = Popen('cat %s/stats_coverage/_%s_all_samples_union_header.bg %s/stats_coverage/%s.bg > %s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir,config_info["abbv"],sp_dir,sample_names_bedgraph_file_list[0],sp_dir,config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
            stdout,stderr = proc.communicate()
            proc = Popen('gzip %s/stats_coverage/_%s_all_samples_union.bg'%(sp_dir,config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
            stdout,stderr = proc.communicate()
            
            #Delete sample bedgraph and header file
            proc = Popen('rm %s/stats_coverage/%s.bg'%(sp_dir,sample_names_bedgraph_file_list[0]),shell=True)
            proc = Popen('rm %s/stats_coverage/_%s_all_samples_union_header.bg'%(sp_dir,config_info["abbv"]),shell=True)



    #Create and submit job to create bedgraph from sum across all sample coverages, compute median, then create clean sites bedgraph, too high sites bedgraph, and too low sites bedgraph, but only if merged bed files don't already exist.

    if os.path.isfile('%s/stats_coverage/_%s_clean_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.isfile('%s/stats_coverage/_%s_too_high_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.isfile('%s/stats_coverage/_%s_too_low_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.getsize('%s/stats_coverage/_%s_clean_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) > 10 and os.path.isfile('%s/stats_coverage/summed_job_complete_no_errors'%sp_dir):
        pass
        
    else:
        sum_coverage_sbatch_file = sum_coverage_sbatch(sp_dir,config_info["abbv"])
        sum_job_id = sbatch_submit(sum_coverage_sbatch_file,memory=8,timelimit=72)
        sleep(30)

        if sum_job_id is not None:
            dones = ['COMPLETED','CANCELLED','FAILED','TIMEOUT','PREEMPTED','NODE_FAIL']
            #Check job id status of sum_job_id job. If not in one of the 'done' job status categories, wait 30 seconds and check again.
            while jobid_status(sum_job_id,start_date) not in dones:
                sleep(30)

            #Check to make sure job completed, and that all necessary files are present. If not, exit and give information.
            sum_job_completion_status = jobid_status(sum_job_id,start_date)
            if sum_job_completion_status != 'COMPLETED':
                sys.exit("There was a problem creating the summed coverage bedgraph file. The job exited with status %s. Please diagnose and fix before moving on"%sum_job_completion_status)
            
            else:
                if os.path.getsize('%s/logs/sumcov_%s.err'%(sp_dir,sum_job_id)) < 1:
                    done_file = open('%s/stats_coverage/summed_job_complete_no_errors'%(sp_dir),"w")
                    done_file.close()
                else:
                    sys.exit('Check whether %s/logs/sumcov_%s.err is empty'%(sp_dir,sum_job_id))



    #Remove unnecessary files if all merged files exist and clean merged file has a filesize > 10 bytes
    if os.path.isfile('%s/stats_coverage/_%s_clean_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.isfile('%s/stats_coverage/_%s_too_high_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.isfile('%s/stats_coverage/_%s_too_low_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) and os.path.getsize('%s/stats_coverage/_%s_clean_coverage_sites_merged.bed'%(sp_dir, config_info["abbv"])) > 10:
        print("All merged bed files exist")
        proc = Popen('rm %s/stats_coverage/_%s_clean_coverage_sites.bed'%(sp_dir, config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
        stdout,stderr = proc.communicate()
        proc = Popen('rm %s/stats_coverage/_%s_too_high_coverage_sites.bed'%(sp_dir, config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
        stdout,stderr = proc.communicate()
        proc = Popen('rm %s/stats_coverage/_%s_too_low_coverage_sites.bed'%(sp_dir, config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
        stdout,stderr = proc.communicate()
        proc = Popen('rm %s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir, config_info["abbv"]),shell=True,stdout=PIPE,stderr=PIPE)
        stdout,stderr = proc.communicate()
    else:
        print('Did not pass all tests to delete files, check if merged bed files exist')





    now = datetime.datetime.now()
    print('Finished script 05: %s'%now)


if __name__ == "__main__":
    main()

