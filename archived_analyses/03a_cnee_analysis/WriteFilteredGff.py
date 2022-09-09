try:
    set
except NameError:
    from sets import Set as set

def ParseMrnaFeatures(annotfeatures):
    rna = annotfeatures.split(';')[0].replace('ID=','') 
    gene = annotfeatures.split(';')[1].replace('Parent=gene-','')
    return rna,gene

def ParseGeneFeatures(annotfeatures):
    gene = annotfeatures.split(';')[0].replace('ID=gene-','')
    return gene

def ParseExonCdsFeatures(annotfeatures):
    #gene = annotfeatures.split(';gene=')[1].split(';')[0]
    mrna = annotfeatures.split(';Parent=')[1].split(';')[0]
    return mrna

pfmrna = open('passfilt_mrna.txt','r')
mrnas = set()
regions = set()
genes = set()
fields = pfmrna.readline().strip().split('\t')
for line in pfmrna:
    linedict = dict(zip(fields,line.strip().split('\t')))
    mrnas.add(linedict['mrnaid'])
    regions.add(linedict['region'])
    genes.add(linedict['geneid'])
gffin = open('galGal6.filt.gff','r')
filtgff = open('galGal6.filtpy.gff','w')
for line in gffin:
    if line[0] == '#':
        filtgff.write(line)
    else:
        linelist = line.strip().split('\t')
        annotfeatures = linelist[8]
        if linelist[2] in ['mRNA','transcript']:
            rna,gene = ParseMrnaFeatures(annotfeatures)
            if rna in mrnas:
                filtgff.write(line)
              
        elif linelist[2] == 'gene':
            gene = ParseGeneFeatures(annotfeatures)
            if gene in genes:
                filtgff.write(line)
        
        elif linelist[2] in ['CDS','exon']:
            mrna = ParseExonCdsFeatures(annotfeatures)
            if mrna in mrnas:
                filtgff.write(line)
        elif linelist[2] == 'region':
            if linelist[0] in regions:
                filtgff.write(line)
        else:
            pass
filtgff.close()

