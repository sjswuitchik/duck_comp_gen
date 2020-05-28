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


#Need:
#sp_dir
#sp_abbr - make abbreviation


#Takes in one line from a union coverage file, and returns a list with the [chromosome,start,end,interval_length,summed_coverage].
def compute_coverage_sum(union_cov_line):
    union_cov_line = union_cov_line.strip()
    union_cov_list = union_cov_line.split()
    
    #Grab positional info, calculate length of interval
    chrom = union_cov_list[0]
    start = union_cov_list[1]
    end = union_cov_list[2]
    interval_length = int(end) - int(start)
    
    #Turn coverage into floats (integers gave problems with e6, etc), sum
    cov_list = [float(i) for i in union_cov_list[3:]]
    summed_coverage = sum(cov_list)
    
    return([chrom,start,end,interval_length,summed_coverage])



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sp_dir", help="directory for that species")
    parser.add_argument("--sp_abbr", help="species abbreviation", required=False)

    args = parser.parse_args()
    sp_dir = args.sp_dir
    sp_abbr = args.sp_abbr
    
    #Read through resulting bedgraph file, create new bedgraph from sum across all sample coverages, and compute median
    #Dictionary to create histogram as we iterate each line
    coverage_histogram = {}    
    total_sites = 0
    if os.path.isfile('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,sp_abbr)) and os.path.getsize('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,sp_abbr)) > 0:
        print("Summary bedgraph and histogram already exists, will not recreate")
   
    else:    

        union_cov_filename = '%s/stats_coverage/_%s_all_samples_union.bg.gz'%(sp_dir, sp_abbr)
        summary_bedgraph = open('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,sp_abbr), 'w')
    
        union_cov_file = gzip.open(union_cov_filename,'r')
        for line in union_cov_file:
            line = line.decode()
            if line[0:5] != 'chrom':
                summed_cov = compute_coverage_sum(line)
                #Returns [chrom,start,end,interval_length,summed_coverage]

                #Write to new summary bedgraph
                summary_bedgraph.write('%s\t%s\t%s\t%d\n'%(summed_cov[0],summed_cov[1],summed_cov[2],summed_cov[4]))
                #Add coverage to histogram dictionary
                if summed_cov[4] not in coverage_histogram:
                    coverage_histogram[summed_cov[4]] = summed_cov[3]
                else:
                    coverage_histogram[summed_cov[4]] = coverage_histogram[summed_cov[4]] + summed_cov[3]
        
                #Add interval length to total sites if coverage is > 0
                if summed_cov[4] > 0:
                    total_sites += summed_cov[3]
 
        union_cov_file.close()
        summary_bedgraph.close()
    
    #If coverage histogram dictionary is empty, open coverage histogram file and read to dictionary, also capture total sites
    if len(coverage_histogram) == 0:
        summed_hist_file = open('%s/stats_coverage/_%s_all_samples_summed_cov_histogram.txt'%(sp_dir,sp_abbr),'r')
        for line in summed_hist_file:
            line = line.strip()
            split_line = line.split()
            if split_line[0] != "SUMMED_COVERAGE":
                coverage_histogram[int(split_line[0])] = int(split_line[1])
                if int(split_line[0]) != 0:
                    total_sites += int(split_line[1])
    
    else:
        #Order and write histogram to file
        ordered_hist_bins = sorted(coverage_histogram.keys())
        summed_hist_file = open('%s/stats_coverage/_%s_all_samples_summed_cov_histogram.txt'%(sp_dir,sp_abbr),'w')
        summed_hist_file.write('SUMMED_COVERAGE\tN_SITES\n')
        for bin in ordered_hist_bins:
            summed_hist_file.write('%d\t%d\n'%(bin,coverage_histogram[bin]))        
            

    #Calculate median value, calculate what the median value is in number of sites, and iterate through bins until that is passed
    ordered_hist_bins = sorted(coverage_histogram.keys())

    median_cov = None
    median_cov_cutoff = (total_sites+1)/2
    cumulative_sites = 0    
    
    for bin in ordered_hist_bins:
        if cumulative_sites < median_cov_cutoff and bin != 0:
            cumulative_sites += coverage_histogram[bin]
            median_cov = bin
        else:
            pass
    
    median_info_file = open('%s/stats_coverage/_%s_median_coverage_info.txt'%(sp_dir,sp_abbr),'w')
    #Calculate upper and lower limits of coverage and print to log file
    upper_lim_cov = median_cov * 2
    lower_lim_cov = median_cov * 0.5
    median_info_file.write('\nMedian coverage of all samples is %d, will exclude all sites > %d and < %d'%(median_cov,upper_lim_cov,lower_lim_cov))
    median_info_file.close()

    #Create new bedfile of sites < 2X median and > 0.5X median
    summary_bedgraph = open('%s/stats_coverage/_%s_all_samples_summed_cov.bg'%(sp_dir,sp_abbr), 'r')
    clean_cov_bedfile = open('%s/stats_coverage/_%s_clean_coverage_sites.bed'%(sp_dir,sp_abbr),'w')
    low_cov_bedfile = open('%s/stats_coverage/_%s_too_low_coverage_sites.bed'%(sp_dir,sp_abbr),'w')
    high_cov_bedfile = open('%s/stats_coverage/_%s_too_high_coverage_sites.bed'%(sp_dir,sp_abbr),'w')
    problem_bedfile = open('%s/stats_coverage/_%s_error_coverage_sites.bed'%(sp_dir,sp_abbr),'w')
    for line in summary_bedgraph:
        line = line.strip()
        split_line = line.split()
        if float(split_line[3]) >= lower_lim_cov and float(split_line[3]) <= upper_lim_cov:
            clean_cov_bedfile.write('%s\t%s\t%s\n'%(split_line[0],split_line[1],split_line[2]))
        elif float(split_line[3]) < lower_lim_cov:
            low_cov_bedfile.write('%s\t%s\t%s\t%s\n'%(split_line[0],split_line[1],split_line[2],split_line[3]))
        elif float(split_line[3]) > upper_lim_cov:
            high_cov_bedfile.write('%s\t%s\t%s\t%s\n'%(split_line[0],split_line[1],split_line[2],split_line[3]))
        else:
            problem_bedfile.write('%s\t%s\t%s\t%s\n'%(split_line[0],split_line[1],split_line[2],split_line[3]))


    summary_bedgraph.close()
    summed_hist_file.close()
    clean_cov_bedfile.close()
    low_cov_bedfile.close()
    high_cov_bedfile.close()
    problem_bedfile.close()
    
if __name__ == "__main__":
    main()

