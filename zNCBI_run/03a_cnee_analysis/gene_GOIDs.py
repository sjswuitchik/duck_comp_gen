#!/usr/bin/env python3

# 7/28/2020 
# ticket INC12977159
# Jacob Pessin    jpessin@bu.edu    Research Computing Services


import sys
from collections import defaultdict


def readfile(fname):
    """Read cness_genes_terms to double-dictionary.

    {result_number: {
        GOID: <goid str value>,
        TERM: <term str value>,
         . . .
        ANNOTATED_GENES: <annotations string> """
    read_results = defaultdict(dict)
    num = None
    with open(fname) as ifile:
        contents = ifile.readlines()

    for line in contents:
        if line.startswith("--"):
            num = line.split()[1]
        elif line.startswith('\n'):
            continue
        elif line and num:  # line is not empty AND past header (first result sub-header)
            line = line.strip("\n")
            key, val = line.split('\t', maxsplit=1) # key (split on first tab)
            if key == 'ANNOTATED_GENES':
                val = val.replace(",", "") # replace comma with nothing
                val = val.split()  # split on whitespace into list
            read_results[num][key] = val
    return read_results


def _flip_goid_gns(doubledict):
    """Created a dict that flips the GOID - GeneName relationship
        --> { Keys: { 'ANNOTATED_GENES': <GENE>, 'GOIDS': [GOID1, 'GOID2', ...] } }
    Note: GOID list is per occurance """
    flipped = defaultdict(list)
    for key in doubledict.keys():
        for gene in doubledict[key]['ANNOTATED_GENES']:
            flipped[gene].append(doubledict[key]['GOID'])
    return flipped


def print_stdout_gn_goids(doubledict):
    """Print to stdout:
        GENE\tGOID;GOID;GOID;... with unique GOIDs"""
    flipped = _flip_goid_gns(doubledict)
    for gene, goids in flipped.items():
        goids = sorted(set(goids)) # set will make unique, sorted returns a list
        goid_str = ";".join(goids)

        line = "\t".join((gene, goid_str))
        print(line)




if __name__ == "__main__":
    contents = readfile(sys.argv[1])
    print_stdout_gn_goids(contents)
