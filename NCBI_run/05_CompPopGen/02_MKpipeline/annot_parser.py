#!/usr/bin/python -tt


# To run:

#   python annot_parser.py $INSHORT.ann.vcf $INSHORT.ann.bed -key missense_variant -key synonymous_variant



import sys
from cyvcf2 import VCF
# see https://brentp.github.io/cyvcf2/docstrings.html

import argparse

#install with "conda install tqdm"
# gives a nice progress bar
import tqdm


def extract_effect(ann, key_effects):
    ''' ann: single input annotation string
            example:  C|synonymous&stop_retained_variant|LOW|id6|GENE_id6|transcript|rna2|....
                  Or:  C|synonymous|LOW|id6|GENE_id6|....
                  Or:  C|stop_retained_variant&synonymous|LOW|.....
        key_effects: list of effect keywords to potentially extract
        output: synonymous 
        
        check out the doctest library in Python
        >>> extract_effect('C|synonymous&stop_retained_variant|LOW|id6|GE',['synonymous', 'missense'])
        ['synonymous']
        
    '''
    # Pull the 2nd element and split by &
    effect_lst = ann.split('|')[1].split('&')
    # check if any key_effects show up in the effect list
    key_set =set(key_effects)
    # do it via set intersection
    effect_set = key_set & set(effect_lst)
    
    # Error check: at most 1 of the key_effects should be in the
    # resulting effect_lst.  If >1, throw an exception
    if len(key_set & effect_set) > 1:
        raise Exception("More than one of the key effects found in the same annotation!")
    
    # return the found effects as a list
    return list(effect_set)



def main():
        parser = argparse.ArgumentParser(description='Do something with VCF files')
        parser.add_argument('infile', help='Input VCF file.')
        parser.add_argument('outfile', help='Output .bed file')
        parser.add_argument('-key',help='Effects to extract. Can be repeated, must be mutually exclusive', action='append')
        args = parser.parse_args()
    
        vcf = VCF(args.infile, gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles

        with open(args.outfile, 'w') as outFile:
            for variant in tqdm.tqdm(vcf):                    
                if len(variant.ALT) == 1:
                    try:
                        for field in variant.INFO:
                                if field[0] == 'ANN':
                                    f1 = field[1]
                                    effect = {}
                                    for ann in f1.split(','):
                                        tmp = extract_effect(ann, args.key)
                                        if tmp:
                                            # Found some!
                                            for found_effect in tmp:
                                                effect[found_effect] = 1
                                    if effect:
                                        # Join the effects after unique-ifying them
                                        if "missense_variant" in effect:
                                            effect_str = "missense_variant"
                                        elif "synonymous_variant" in effect:
                                            effect_str = "synonymous_variant"
                                        else:
                                            raise Exception("Found an effect but it doesn't contain any keys!")
                                        start_pos = int(variant.end) - 1
                                        outFile.write('%s\t%s\t%s\t%s\n' % (variant.CHROM, start_pos, variant.end, effect_str))
                    except Exception as e:
                        # two keys found in the same annotation.
                        sys.stderr('Two mutually exclusive effects found at the same time.')
                        sys.stderr(variant)
                        
            outFile.close()
            
            
    
if __name__ == '__main__':
  main()


