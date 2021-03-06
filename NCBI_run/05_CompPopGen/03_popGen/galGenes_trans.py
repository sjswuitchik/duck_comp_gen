import fileinput
import sys

def gene_abrvs(line):
    """Abbreviations imbedded as ;gene=CONTENT"""
    parts = line.split(';')
    abrvs = [part.rsplit('=')[1] for part in parts if 'gene=' in part]
    abrvs = set(abrvs) # uniq, arbitrary order
    return abrvs
    
    
def id_abrvs(line):
    """Abbreviations imbedded as ;ID=rna"""
    parts = line.split(';')
    numabrvs = [part.rsplit('=')[1] for part in parts if 'ID=rna' in part]
    numabrvs = set(numabrvs)
    return numabrvs


with fileinput.input() as intake:
    for line in intake:
        abrvset = gene_abrvs(line)
        abrvset = ",".join(list(set(abrvset)))
        numabrvset = id_abrvs(line)
        numabrvset = ",".join(list(set(numabrvset)))
        result = '\t'.join((abrvset, numabrvset))
        print(result, end='\n', file=sys.stdout)
