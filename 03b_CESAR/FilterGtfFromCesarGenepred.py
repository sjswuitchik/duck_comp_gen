## filtering a sorted GTF created by GenePredToGtf from a CESAR 2.0 GenePred ## 
# AH Freedman, May 2020 # 

import argparse

def BuildScaffoldSizeDictionary(infile):
    size_dict = {}
    fopen = open(infile,'r')
    for line in fopen:
        id,size=line.strip().split('\t')
        size_dict[id] = int(size)
    return size_dict

def ParseGtfFeature(featureline,gtf_dict):
    scaffold,source,feature,start,end,score,strand,frame,attribute = featureline.strip().split('\t')
    if feature not in ['CDS','exon','start_codon','transcript','stop_codon']:
        raise Exception('unrecognized feature class')
    else:
        attribute_dict = {}
        for entry in attribute.replace('"','').split(';')[:-1]:
            attribute_dict[entry.split()[0]] = entry.split()[1]

        if attribute_dict['transcript_id'] not in gtf_dict:
            gtf_dict[attribute_dict['transcript_id']] = {}
            gtf_dict[attribute_dict['transcript_id']]['scaffolds'] = []
            gtf_dict[attribute_dict['transcript_id']]['lines'] = []
            gtf_dict[attribute_dict['transcript_id']]['starts'] = []
            gtf_dict[attribute_dict['transcript_id']]['ends'] = []
            gtf_dict[attribute_dict['transcript_id']]['features'] = []
                
        gtf_dict[attribute_dict['transcript_id']]['scaffolds'].append(scaffold)
        gtf_dict[attribute_dict['transcript_id']]['lines'].append(featureline.strip())
        gtf_dict[attribute_dict['transcript_id']]['starts'].append(int(start))
        gtf_dict[attribute_dict['transcript_id']]['ends'].append(int(end))
        gtf_dict[attribute_dict['transcript_id']]['features'].append(feature)
    return gtf_dict,attribute_dict['transcript_id']
    
def RemoveBadIntervals(gtf_dict,scaffold_size_dict):
    filtered_dict = {}
    for key in gtf_dict:
        flag = 0
        for i in range(len(gtf_dict[key]['features'])):
            if gtf_dict[key]['starts'][i] >= gtf_dict[key]['ends'][i]:
                flag+=1
            elif gtf_dict[key]['starts'][i] > scaffold_size_dict[gtf_dict[key]['scaffolds'][i]]:
                flag+=1
            elif gtf_dict[key]['ends'][i] > scaffold_size_dict[gtf_dict[key]['scaffolds'][i]]:
                flag+=1
        if flag == 0:
            filtered_dict[key] = gtf_dict[key]
        
    return filtered_dict

def RestoreAnnotationHierarchy(gtf_dict,hierarchy=None):
    hierarchy = ['transcript','exon','CDS','start_codon','stop_codon']
    for transcript in gtf_dict:
        reordered_features = []
        for level in hierarchy:
            for i in range(len(gtf_dict[transcript]['features'])):
                if gtf_dict[transcript]['features'][i] == level:
                    reordered_features.append(gtf_dict[transcript]['lines'][i])
        gtf_dict[transcript]['reordered_features'] = reordered_features

    return gtf_dict


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='options for filtering gtf derived from CESAR2.0 genepred-format annotation')
    parser.add_argument('-gtf','--cesar-gtf',dest='gtf',type=str,help='sorted gtf converted from CESAR2.0 genepred annotation')
    parser.add_argument('-sizes','--genome-scaffold-sizes',dest='sizes',type=str,help='table of scaffold lengths')
    opts = parser.parse_args()

    scaffold_size_dict = BuildScaffoldSizeDictionary(opts.sizes)
    gtfopen = open(opts.gtf,'r')
    tscript_ids = []
    gtf_dict = {}
    for line in gtfopen:
        gtf_dict,transcript_id = ParseGtfFeature(line,gtf_dict)
        if transcript_id not in tscript_ids:
            tscript_ids.append(transcript_id)

    interval_cleaned_gtf_dict = RemoveBadIntervals(gtf_dict,scaffold_size_dict)
    reordered_cleaned_gtf_dict = RestoreAnnotationHierarchy(interval_cleaned_gtf_dict)
    reordered_cleaned_out = open('cleaned_reordered_%s' % opts.gtf,'w')
    for transcript in tscript_ids:
        if transcript in reordered_cleaned_gtf_dict:
            reordered_cleaned_out.write('%s\n' % '\n'.join(reordered_cleaned_gtf_dict[transcript]['reordered_features']))

reordered_cleaned_out.close()
