from collections import defaultdict
try:
    set
except NameError:
    from sets import Set as set

def ExtractMrnaId(iddata):
    rna = iddata.split(';')[0].replace('ID=','') 
    gene = iddata.split(';')[1].replace('Parent=gene-','')
    return rna,gene

gffin = open('galGal6.filt.gff','r')

mrnaout = open('passfilt_mrna.txt','w')
mrnaout.write('region\tmrnaid\tgeneid\n')
genes = set()
for line in gffin:
    if line[0] != '#':
        linelist = line.strip().split('\t')
        if linelist[2] == 'gene' and 'gene_biotype=protein_coding' in line:
            genes.add(linelist[8].split(';')[0].replace('ID=gene-','')) 
        elif linelist[2] == 'mRNA' and 'partial=true' not in line or linelist[2] == 'transcript':
            annotdata = linelist[8]
            rna,gene=ExtractMrnaId(annotdata)
            if gene in genes:
                if 'frameshifting' not in annotdata:
                    mrnaout.write('%s\t%s\t%s\n' % (linelist[0],rna,gene))
                elif 'non-frameshifting' in annotdata and ' frameshifting' not in annotdata:
                    mrnaout.write('%s\t%s\t%s\n' % (linelist[0],rna,gene))


mrnaout.close()

