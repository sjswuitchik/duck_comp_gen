#!/usr/bin/env python3

import os, sys, argparse

print()

#set up parser for collecting input parameters through command line
parser = argparse.ArgumentParser(description=
                                 'This converts a newick tree file to a nexus tree file.')
requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type=str, metavar='input_newick', required=True, help='Name of input newick file')
requiredParam.add_argument('-o', type=str, metavar='output_nexus', required=True, help='Name of output nexus file')
args = parser.parse_args()

infile = open(args.i,'r')
outfile = open(args.o,'w')
outfile.write('#NEXUS\n\nBegin trees;\n')

count = 0
tree = ''

print('Converting trees...\n')
for line in infile:
    data = line.strip('\n')
    tree = tree+data
    if data.endswith(';'):
        count += 1
        outfile.write('tree tree'+str(count)+' = [&U] '+tree+'\n')
        tree = ''
outfile.write('end;')

infile.close()
outfile.close()

print(str(count)+' trees converted.\n\nFinished!!\n\n')

