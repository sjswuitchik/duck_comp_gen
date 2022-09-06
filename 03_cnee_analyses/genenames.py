import fileinput
import sys

def gene_abrvs(line):
    """Abbreviations imbedded as ;gene=CONTENT"""
    parts = line.split(';')
    abrvs = [part.rsplit('=')[1] for part in parts if 'gene=' in part]
    abrvs = set(abrvs) # uniq, arbitrary order
    return abrvs


def location(line):
    """Extracts chrom, chromStart, chromStop (columns 1-3 in a .bed file)
		returns as a tab seperated string"""
    poss = line.split()[0:3]
    poss = "\t".join(poss)
    return poss


with fileinput.input() as intake:
    for line in intake:
        abrvset = gene_abrvs(line)
        abrvset = ",".join(list(set(abrvset)))
        loc = location(line)
        result = '\t'.join((loc, abrvset))
        print(result, end='\n', file=sys.stdout)
