#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import shutil
import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg') #Added to get working on the cluster
from matplotlib.backends.backend_pdf import PdfPages


#####May wish to run scripts 01 and 02 on regal until final deduplicated BAM is produced. After that, copy to "permanent" directory on holylfs. Create script 2.5 to just set up working dir on holylfs, mv dedup dir, logs, scripts, stats over and end.

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
    
    
##########Above functions same as script 01, may just wish to source that script in future.


#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def dedup_sbatch(sp_dir,sp_abbr,sample_ncbi_dict):
    slurm_script = script_create()
    dedup_sbatch_filenames = []
     
    for sample in sample_ncbi_dict.keys():
        #First check if dedup file is already present (already downloaded), or final BAM file already present. If it has, print statment and continue with next sample. 
        dedup_filename = '%s/dedup/%s.dedup.bam'%(sp_dir,sample)
        dedup_sorted_filename = '%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)
        if os.path.isfile(dedup_filename) or os.path.isfile(dedup_sorted_filename):
            print('%s.dedup.bam already present, skipping'%(sample))
        else:
            print('Will dedup sras for sample %s'%(sample))
            

            align_dir = '%s/alignment/'%(sp_dir)

            #Load modules and get versions for all programs used
            ##For now, using my own installation of GATK as it is not yet installed on the cluster
            cmd_1 = 'module load jdk/1.8.0_45-fasrc01'
            

            #Add PL read group info (can remove once all initial libraries done with script 2)
            rg_cmd_list = []
            for sra in sample_ncbi_dict[sample]:
                
                rg_cmd = '/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05_CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" AddOrReplaceReadGroups -I %s/alignment/%s.sorted.bam -O %s/alignment/%s.sorted.rg.bam -ID %s -SM %s -PU %s.%s -LB %s -PL illumina --COMPRESSION_LEVEL 5 --CREATE_INDEX true'%(sp_dir,sra,sp_dir,sra,sra,sample,sra,sample,sample)
                
                rg_cmd_list.append(rg_cmd)
            
            rg_cmd = "\n\n".join(rg_cmd_list)
            
            #Create GATK mark duplicates command
            gatk_cmd_1 = '/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05_CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" MarkDuplicatesGATK '
            gatk_cmd_2 = ('-I '+align_dir+'%s.sorted.rg.bam ')*len(sample_ncbi_dict[sample])%tuple(sample_ncbi_dict[sample])
            gatk_cmd_3 = '-O %s/dedup/%s.dedup.bam '%(sp_dir,sample)
            gatk_cmd_4 = '--METRICS_FILE %s/stats/%s.dedup.metrics.txt '%(sp_dir,sample)
            gatk_cmd_5 = '--COMPRESSION_LEVEL 5'
            
            #Combine into single line
            cmd_2 = "".join([gatk_cmd_1,gatk_cmd_2,gatk_cmd_3,gatk_cmd_4,gatk_cmd_5])
            
            #Sort and index dedup bam file
            cmd_3 = '/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05_CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" SortSam -I %s/dedup/%s.dedup.bam -O %s/dedup/%s.dedup.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --COMPRESSION_LEVEL 5'%(sp_dir,sample,sp_dir,sample)
            
            cmd_4 = '/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05_CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" CollectAlignmentSummaryMetrics -I %s/dedup/%s.dedup.sorted.bam -R %s/genome/%s.fa --METRIC_ACCUMULATION_LEVEL=SAMPLE -O %s/stats/%s.alignment_metrics.txt'%(sp_dir,sample,sp_dir,sp_abbr,sp_dir,sample)
            
            #Validate sorted bam
            cmd_5 = '/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/05_CompPopGen/gatk-4.0.3.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=1" ValidateSamFile -I %s/dedup/%s.dedup.sorted.bam -O %s/stats/%s.validate.txt'%(sp_dir,sample,sp_dir,sample)
            
            #Compute coverage histogram of sorted bam
            cmd_6 = 'bedtools genomecov -ibam %s/dedup/%s.dedup.sorted.bam -g %s/genome/%s.fa > %s/stats/%s.coverage'%(sp_dir,sample,sp_dir,sp_abbr,sp_abbr,sample)
            
            #Grab only genome output
            cmd_7 = r"""awk '$1 == "genome" {print $0}' %s/stats/%s.coverage > %s/stats/%s.genome.coverage"""%(sp_dir,sample,sp_dir,sample)
        
            cmd_list = [cmd_1,rg_cmd,cmd_2,cmd_3,cmd_4,cmd_5,cmd_6,cmd_7]

            final_cmd = "\n\n".join(cmd_list)


    #Format sbatch script
            sra_script = slurm_script.format(partition="shared",time="1-12:00",mem="10000",cores="2",nodes="1",jobid="dedup",sp_dir=sp_dir,cmd=final_cmd)
            out_filename = "%s/scripts/04_dedup_sort_validate_%s.sbatch"%(sp_dir,sample)
            out_file = open(out_filename,"w")
            out_file.write(sra_script)
            out_file.close
            dedup_sbatch_filenames.append(out_filename)
    
    return(dedup_sbatch_filenames)


#Collect relevant line of alignment stats (from pair) for each sample. Returns a dictionary of header column names as keys with values for that sample as values.
def collect_alignment_metrics(sp_dir,sample):
    align_stats_file = '%s/stats/%s.alignment_metrics.txt'%(sp_dir,sample)
    align_stats = {}
    try:
        full_align_stats = open(align_stats_file,"r")
        for line in full_align_stats:
            split_line = line.strip().split("\t")
            if split_line[0] == "CATEGORY":
                header = split_line
            #Look for lines with the paired data at the sample level (no read group but sample in sample list, return relevant line)
            else:
                try:
                    if split_line[24] == sample and split_line[0] == "PAIR":
                        align_stats = dict(zip(header,split_line))
                except:
                    pass
    except:
        print('No alignment stats file for sample: %s'%(sample))

    return(align_stats)
    
#Collect all deduplication metrics for a list of sample IDs, writes all metrics to a file and returns a dictionary of % duplication for each sample
def collect_dedup_metrics(sp_dir,sample):
    dedup_stats_file = '%s/stats/%s.dedup.metrics.txt'%(sp_dir,sample)
    dedup_stats = {}
    try:
        full_dedup_stats = open(dedup_stats_file,"r")
        for line in full_dedup_stats:
            split_line = line.strip().split("\t")
            if split_line[0] == "LIBRARY":
                header = split_line
            #Look for lines with the paired data at the sample level (no read group but sample in sample list, return relevant line)
            elif split_line[0] == sample:
                dedup_stats = dict(zip(header,split_line))
    except:
        print('No dedup stats file for sample: %s'%(sample))

    return(dedup_stats)
    
#Given a sample name, reads in histogram of genome coverage, calculates mean and median, and returns dictionary with coverage histogram, mean and median for that sample
def collect_coverage_metrics(sp_dir,sample):
    coverage_file = '%s/stats/%s.genome.coverage'%(sp_dir,sample)
    coverage_stats = {}
    hist_bins = []
    hist_vals = []
    site_total = 0
    cov_total = 0
    try:
        full_coverage_stats = open(coverage_file,"r")
        for line in full_coverage_stats:
            split_line = line.strip().split("\t")
            med_cutoff = (int(split_line[3])+1)/2
            site_total += int(split_line[2])
            cov_total += (int(split_line[1])*int(split_line[2]))
            
            hist_bins.append(int(split_line[1]))
            hist_vals.append(float(split_line[4]))
            
            #Test whether site total is greater than midpoint if median not yet discovered - if yes, add to median in dictionary.
            if "median" not in coverage_stats and site_total > med_cutoff:
                coverage_stats["median"] = int(split_line[1])
            else:
                pass
        
        #Calculate mean from cov_total and site_total
        coverage_stats["mean"] = round(cov_total/site_total,2)
        
        #Add histogram to dictionary
        coverage_stats["hist_vals"] = hist_vals
        coverage_stats["hist_bins"] = hist_bins
    except:
        print("No genome coverage file for sample: %s"%(sample))
        
    return(coverage_stats)
    
    


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping")
    args = parser.parse_args()
    config_filename = args.config
    
    now = datetime.datetime.now()
    print('Staring work on script 02: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get Sample, SRA and Genome attributes - use same config as for pipeline script 01    
    config_info = extract_config(config_filename)


    #####Check if dedup directory exists, if not creates it
    
    sp_dir = "%s/%s"%(config_info["out_dir"],config_info["abbv"])
    
    dedup_dir = "%s/dedup"%(sp_dir)
    directory_create(dedup_dir)


    #####Create sbatch files to dedup SRA files and combine if multiple SRAs for each sample, sort and index resulting file, validate, and compute coverage histogram
    
    #Create sbatch files
    dedup_filenames = dedup_sbatch(sp_dir,config_info["abbv"],config_info["sample_dict"])
    
    #Submit dedup read sbatch files
    dedup_jobids = []
    completed_jobids = {}
    for i in range(0,len(dedup_filenames)):
        dedup_jobids.append(sbatch_submit(dedup_filenames[i]))
        sleep(1)
    #Add an extra sleep to give sacct a chance to catch up
    sleep(20)
    #Then, enter while loop that will continue until the number of completed jobs matches the number of sbatch files
    while len(completed_jobids) < len(dedup_filenames):
        num_running = num_pend_run(dedup_jobids,start_date)
        job_statuses = all_jobs_status(start_date)
        for job in dedup_jobids:
            if job not in completed_jobids:
                if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                    completed_jobids[job] = job_statuses[job]
                    print("Job %s completed"%job)
        sleep(30)
    
    #After all jobs have finished, report which jobs failed
    for job in completed_jobids:
        if completed_jobids[job] != "COMPLETED":
            print("Dedup, sort and validate job %s failed with code: %s"%(job,completed_jobids[job]))
 
    
    #####Collate summary statistics   
    
    #Get stats directory filenames
    stat_files = os.listdir("%s/stats"%(sp_dir))
    
    #Collect alignment metrics into a single file for all samples - if none exists will print statement and add empty list to dictionary
    all_align_stats = {}
    for sample in config_info["sample_dict"]:
        all_align_stats[sample] = collect_alignment_metrics(sp_dir,sample)
        
    #Write all sample alignment metrics to a file
    align_file = open('%s/stats/_%s_all_sample_alignment_stats.txt'%(sp_dir,config_info["abbv"]),"w")
    align_file.write("SAMPLE\tCATEGORY\tTOTAL_READS\tPF_READS\tPCT_PF_READS\tPF_NOISE_READS\tPF_READS_ALIGNED\tPCT_PF_READS_ALIGNED\tPF_ALIGNED_BASES\tPF_HQ_ALIGNED_READS\tPF_HQ_ALIGNED_BASES\tPF_HQ_ALIGNED_Q20_BASES\tPF_HQ_MEDIAN_MISMATCHES\tPF_MISMATCH_RATE\tPF_HQ_ERROR_RATE\tPF_INDEL_RATE\tMEAN_READ_LENGTH\tREADS_ALIGNED_IN_PAIRS\tPCT_READS_ALIGNED_IN_PAIRS\tBAD_CYCLES\tSTRAND_BALANCE\tPCT_CHIMERAS\tPCT_ADAPTER\n")
    for sample in sorted(config_info["sample_dict"].keys()):
        if all_align_stats[sample]:
            sample_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(sample,all_align_stats[sample]["CATEGORY"],all_align_stats[sample]["TOTAL_READS"],all_align_stats[sample]["PF_READS"],all_align_stats[sample]["PCT_PF_READS"],all_align_stats[sample]["PF_NOISE_READS"],all_align_stats[sample]["PF_READS_ALIGNED"],all_align_stats[sample]["PCT_PF_READS_ALIGNED"],all_align_stats[sample]["PF_ALIGNED_BASES"],all_align_stats[sample]["PF_HQ_ALIGNED_READS"],all_align_stats[sample]["PF_HQ_ALIGNED_BASES"],all_align_stats[sample]["PF_HQ_ALIGNED_Q20_BASES"],all_align_stats[sample]["PF_HQ_MEDIAN_MISMATCHES"],all_align_stats[sample]["PF_MISMATCH_RATE"],all_align_stats[sample]["PF_HQ_ERROR_RATE"],all_align_stats[sample]["PF_INDEL_RATE"],all_align_stats[sample]["MEAN_READ_LENGTH"],all_align_stats[sample]["READS_ALIGNED_IN_PAIRS"],all_align_stats[sample]["PCT_READS_ALIGNED_IN_PAIRS"],all_align_stats[sample]["BAD_CYCLES"],all_align_stats[sample]["STRAND_BALANCE"],all_align_stats[sample]["PCT_CHIMERAS"],all_align_stats[sample]["PCT_ADAPTER"])
            align_file.write(sample_line)
        else:
            align_file.write('%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n'%(sample))
    align_file.close()

    
    #Collect dedup metrics into a single file for all samples
    all_dedup_stats = {}
    for sample in config_info["sample_dict"]:
        all_dedup_stats[sample] = collect_dedup_metrics(sp_dir,sample)
        
    #Write all sample dedup metrics to a file
    dedup_file = open('%s/stats/_%s_all_sample_dedup_stats.txt'%(sp_dir,config_info["abbv"]),"w")
    dedup_file.write("SAMPLE\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n")
    for sample in sorted(config_info["sample_dict"].keys()):
        if all_dedup_stats[sample]:
            try:
                sample_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(sample,all_dedup_stats[sample]["UNPAIRED_READS_EXAMINED"],all_dedup_stats[sample]["READ_PAIRS_EXAMINED"],all_dedup_stats[sample]["UNMAPPED_READS"],all_dedup_stats[sample]["UNPAIRED_READ_DUPLICATES"],all_dedup_stats[sample]["READ_PAIR_DUPLICATES"],all_dedup_stats[sample]["READ_PAIR_OPTICAL_DUPLICATES"],all_dedup_stats[sample]["PERCENT_DUPLICATION"],all_dedup_stats[sample]["ESTIMATED_LIBRARY_SIZE"])
                dedup_file.write(sample_line)
            except:
                dedup_file.write('%s\t\t\t\t\t\t\t\t\n'%(sample))
        else:
            dedup_file.write('%s\t\t\t\t\t\t\t\t\n'%(sample))
    dedup_file.close()

    
    #Collect coverage histograms, plot in pdf (species_date.pdf), calculate mean and median coverage
    all_coverage_stats = {}
    for sample in config_info["sample_dict"]:
        all_coverage_stats[sample] = collect_coverage_metrics(sp_dir,sample)
        
    #Plot all coverage histograms in pdf, initialize file
    pdf_file = PdfPages('%s/stats/_%s_all_coverage_plots.pdf'%(sp_dir,config_info["abbv"]))
    
    #Iterate through samples and add plots to pdf file
    for sample in config_info["sample_dict"]:
        if "hist_vals" in all_coverage_stats[sample]:
    
            #Change size of figure 
            plt.figure(figsize=(8,4))
            
            #Only plot up to 3x the mean
            bars = []
            height = []
            for i in range(0,len(all_coverage_stats[sample]["hist_bins"])):
                if all_coverage_stats[sample]["hist_bins"][i] < (all_coverage_stats[sample]["mean"]*3):
                    bars.append(all_coverage_stats[sample]["hist_bins"][i])
                    height.append(all_coverage_stats[sample]["hist_vals"][i])
                    
            y_pos = np.arange(len(bars))
            
            #Create bars
            plt.bar(y_pos,height,color="c")
            
            #Add title, axis names
            plt.title("%s"%sample)
            plt.xlabel("coverage")
            plt.ylabel("proportion reads")
            
            #Change ticks to coverage labels, rotate
            plt.xticks(y_pos,bars, rotation=90, size=6)
            
            #Save plot to pdf
            #plt.show()
            pdf_file.savefig()
    
    pdf_file.close()
    
    #Check the output of validation file. Add to dictionary of samples, with "ok" if everything looks okay, and "check" if there are any issues.
    all_validation_stats = {}
    for sample in config_info["sample_dict"]:
        try:
            valid_file = open('%s/stats/%s.validate.txt'%(sp_dir,sample),"r")
            valid_contents = []
            for line in valid_file:
                line = line.strip()
                valid_contents.append(line)
            if valid_contents[0] == "No errors found":
                all_validation_stats[sample] = "ok"
            else:
                all_validation_stats[sample] = "check"
                print("Sample %s has an error in the validation file"%sample)
        except:
            all_validation_stats[sample] = "check"
            print("Sample %s has an error, or is missing the validation file"%sample)
    
    #Produce file of all most important summary stats (including validation status). 
    summary_stat_file = '%s/stats/_%s_all_summary_stats.txt'%(sp_dir,config_info["abbv"])
    sum_stat = open(summary_stat_file,"w")
    
    #Write header
    sum_stat.write("SAMPLE\tMEAN_COVERAGE\tMEDIAN_COVERAGE\tTOTAL_READS\tMEAN_READ_LENGTH\tPCT_PF_READS_ALIGNED\tPCT_PF_HQ_ALIGNED_READS\tPF_HQ_MEDIAN_MISMATCHES\tPF_INDEL_RATE\tPCT_READS_ALIGNED_IN_PAIRS\tSTRAND_BALANCE\tPERCENT_DUPLICATION\tVALIDATION_STATUS\n")
    
    #Iterate through samples and add results from all three dictionaries if present
    #Also add median coverage to list, to add decision about pipeline recommendation for haplotypecaller
    all_median_cov = []
    
    for sample in config_info["sample_dict"]:
        samp_sum = []
        samp_sum.append(sample)
        #coverage
        if "mean" in all_coverage_stats[sample]:
            samp_sum.append(str(all_coverage_stats[sample]["mean"]))
            samp_sum.append(str(all_coverage_stats[sample]["median"]))
            all_median_cov.append(all_coverage_stats[sample]["median"])
        else:
            samp_sum.append("")
            samp_sum.append("")
        #alignment
        if "TOTAL_READS" in all_align_stats[sample]:
            samp_sum.append(all_align_stats[sample]["TOTAL_READS"])
            samp_sum.append(all_align_stats[sample]["MEAN_READ_LENGTH"])
            samp_sum.append(all_align_stats[sample]["PCT_PF_READS_ALIGNED"])
            samp_sum.append(str(round(int(all_align_stats[sample]["PF_HQ_ALIGNED_READS"])/int(all_align_stats[sample]["PF_READS"]),6)))
            samp_sum.append(all_align_stats[sample]["PF_HQ_MEDIAN_MISMATCHES"])
            samp_sum.append(all_align_stats[sample]["PF_INDEL_RATE"])
            samp_sum.append(all_align_stats[sample]["PCT_READS_ALIGNED_IN_PAIRS"])
            samp_sum.append(all_align_stats[sample]["STRAND_BALANCE"])
        else:
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")
            samp_sum.append("")           
        #dedup
        if "PERCENT_DUPLICATION" in all_dedup_stats[sample]:
            samp_sum.append(all_dedup_stats[sample]["PERCENT_DUPLICATION"])
        else:
            samp_sum.append("")
        
        samp_sum.append(all_validation_stats[sample])
            
        samp_sum = "\t".join(samp_sum)
        sum_stat.write('%s\n'%samp_sum)
    
    sum_stat.close()
   
    #Copy this file and pdf of coverage to centralized location        
    general_dir = "_ALL_SPECIES_SUMMARIES"
    directory_create(general_dir)
    
    try:
        proc = Popen('cp %s %s/%s_all_summary_stats.txt'%(summary_stat_file,general_dir,config_info["abbv"]),shell=True)
        proc = Popen('cp %s/stats/_%s_all_coverage_plots.pdf %s/%s_all_coverage_plots.pdf'%(sp_dir,config_info["abbv"],general_dir,config_info["abbv"]),shell=True)
    except:
        print("There was an error copying summary stat files")
    
    #Extract min and max median coverage, if both fall below 15X, recommend lowcoverage pipeline for haplotype caller. If both fall above 15X, recommend highcoverage pipeline for haplotype caller. If they span 15X, recommend that the user examine coverages to determine best pipeline to use.
    
    max_median_cov = max(all_median_cov)
    min_median_cov = min(all_median_cov)
    
    if max_median_cov > 15 and min_median_cov > 15:
        print("All samples median coverage > 15X, recommended pipline for script 03: highcoverage\n")
    elif max_median_cov < 15 and min_median_cov < 15:
        print("All samples median coverage < 15X, recommended pipline for script 03: lowcoverage\n")
    elif max_median_cov > 15 and min_median_cov < 15:
        print("Some samples median coverage > 15X and some < 15X, check summary stats file to determine best pipeline\n")
    
    #Check that the final sorted bam and index is available, if so, remove intermediate files (unsorted dedup, all aligned bams)
    ###Only delete if validate says no errors found.
        ####Move this to script 2 after validation has occurred.         
    #Check that the final sorted bam and index is available, if so, remove intermediate files (SRA, fastq,unsorted BAM)

    for sample in config_info["sample_dict"]:
        if os.path.isfile('%s/dedup/%s.dedup.sorted.bam'%(sp_dir,sample)) and os.path.isfile('%s/dedup/%s.dedup.sorted.bai'%(sp_dir,sample)) and all_validation_stats[sample] == "ok":
            if os.path.isfile('%s/dedup/%s.dedup.bam'%(sp_dir,sample)):
                proc = Popen('rm %s/dedup/%s.dedup.bam'%(sp_dir,sample),shell=True)
            for sra in config_info["sample_dict"][sample]:
                if os.path.isfile('%s/alignment/%s.sorted.rg.bam'%(sp_dir,sra)):
                    proc = Popen('rm %s/alignment/%s.sorted*'%(sp_dir,sra),shell=True)
                if os.path.isfile('%s/fastq/%s_1.fastq.gz'%(sp_dir,sra)):
                    proc = Popen('rm %s/fastq/%s*'%(sp_dir,sra),shell=True)
                if os.path.isfile('%s/sra/%s.sra'%(sp_dir,sra)):
                    proc = Popen('rm %s/sra/%s.sra'%(sp_dir,sra),shell=True)
        else:
            print("Something happened with sample deduping: %s"%(sample))     
               
    now = datetime.datetime.now()
    print('Finished script 02: %s'%now)

if __name__ == "__main__":
    main()

