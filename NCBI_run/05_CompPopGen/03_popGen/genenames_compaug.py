import fileinput
import sys

def gene_abrvs(line):
    """Abbreviations imbedded as gene_id "CONTENT"; """
    parts = line.split()
    abrv = parts[parts.index("gene_id")+1][1:-2] if "gene_id" in parts else "";
    abrv = set([] if abrv == "" else [abrv])
    return abrv

def location(line):
    """Extracts chrom, chromStart, chromStop (columns 1-3 in a .bed file)
		returns as a tab seperated string"""
    poss = line.split()[0,4:5]
    poss = "\t".join(poss)
    return poss


with fileinput.input() as intake:
    for line in intake:
        abrvset = gene_abrvs(line)
        abrvset = ",".join(list(set(abrvset)))
        loc = location(line)
        result = '\t'.join((loc, abrvset))
        print(result, end='\n', file=sys.stdout)
