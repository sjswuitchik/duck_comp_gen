#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
import datetime

#Check if a species directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):
        os.makedirs(dir)

#Extract info from local fastq config file    
def parse_local_config(config_filename):
    #local config has four columns: path to fastq, read (#1 or #2), sample ID, "accession"    
    with open(config_filename) as f:
        lines = f.read().splitlines()
    return(lines)
    
#Write config file
def write_main_config(config_local, config_out, abbv):
    with open(config_out, mode="w") as f:
        print("#Example config file:", file=f)
        print("", file=f)
        print("#Path to desired output directory (directory for each species will be created within this directory). If not specified will be current working directory", file=f)
        print("--OUT_DIR .", file=f)
        print("", file=f)
        print("#Species abbreviation for use in file headers (e.g. CCornix for hooded crow)", file=f)
        print("--ABBV", abbv, sep=" ", file=f)
        print("", file=f)
        
        #now do each sample_local line
        for line in config_local:
            line=line.split()
            #path = line[0], read = line[1], sample ID = line[2], accession = line[3]
            if line[1] == '1':
                print("--SAMPLE_LOCAL", line[2], line[3], sep=" ", file=f)
    
#Copy and rename fastqs
def cp_rename_fastq(config_local, sp_dir):
    #takes the local config dict and the species dir
    #checks if desired fastq.gz files already exist, and if not:
    #generates commands to copy the local fastq to species_dir/ACCESSION_1.fastq.gz
    #checks to make sure that the input is actually gzipped only by checking for extension
    #will gzip fastqs if not already, otherwise, just copies them
    for line in config_local:
        line = line.split()
        #path = line[0], read = line[1], sample ID = line[2], accession = line[3]
        #if path ends with .gz, use gzip -cd
        #else use cp
        if os.path.isfile("%s/fastq/%s_%s.fastq.gz"%(sp_dir,line[3],line[1])) is False:
            if line[0].endswith(".gz"):
                command = "cp %s %s/fastq/%s_%s.fastq.gz"%(line[0],sp_dir,line[3],line[1])
            else:
                command = "gzip -c %s > %s/fastq/%s_%s.fastq.gz"%(line[0],sp_dir,line[3],line[1])
            result=os.system(command)
        else:
            pass


def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples and genome for mapping", required=True)
    parser.add_argument("--out_dir", help="path to the desired output directory", required=True)
    parser.add_argument("--abbv", help="species abbreviation to use", required=True)
    parser.add_argument("--config_out", help="write config file for full pipeline to this file", required=True)
    args = parser.parse_args()
    
    now = datetime.datetime.now()
    print('Staring work: %s'%now)
    start_date = now.strftime("%Y-%m-%d")

    config_local = parse_local_config(args.config)
    sp_dir = "%s/%s"%(args.out_dir,args.abbv)
    print("\nOutput will be written to %s\n"%sp_dir)
    directory_create(sp_dir)
    directory_create("%s/fastq"%(sp_dir))
    #actually run cp
    cp_rename_fastq(config_local, sp_dir)
    #make config for main pipeline
    write_main_config(config_local, args.config_out, args.abbv)
    
if __name__ == "__main__":
    main()

