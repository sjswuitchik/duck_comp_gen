#!/bin/bash
# NAME
#     generate_orthogroups_fasta.sh - generate one FASTA per orthogroup in current working directory
#
# USAGE
#     cd output_dir 
#     /path/to/generate_orthogroups_fasta.sh
#
# INPUT FILES
#     orthogroup TSV file of the form:
#
#     HEADER LINE (unused...)
#     OG0001593	NP_001297327.1, XP_027323025.1, XP_027323027.1	NP_001297327.1, XP_027323025.1	jg19339.t1    ...
#     OG0001594	NP_001297338.1, XP_021132300.2, XP_027316613.1, XP_027317818.1	NP_001297338.1, XP_027316613.1	...
#     ...
#     
#     Currently, this script skips the header line and 2nd (tab-separated) column
#     Columns are assumed to be ordered to correspond to "species" and species-fasta-file arguments to the awk script below...
#    
# OUTPUT FILES
#     One FASTA file per orthogroup; e.g.:
#     $ head OG0001593_nuc.fa
#     >anaPla_NP_001297327.1
#     ATGGTTGACACAGAA...
#     >anaPla_XP_027323025.1
#     ATGGTTGACACAGAA...
#     >ansBra_jg19339.t1
#     ATGGTTGACACAGAA...
#     >ansCyg_XP_013043163.1
#
# STDERR
#     Duplicate sequence IDs detected in the input FASTA files are detected and reported.
#     The first sequence encountered for a given seqid is used in the output.
#
# AUTHOR
#     Nathan Weeks <nweeks@g.harvard.edu>
#
# VERSION HISTORY
#     2021-11-30    Initial version

set -o errexit -o nounset -o pipefail

# remove header
tail -n +2 /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/clean_orthos.tsv |
  # remove (unused) 2nd column
  cut -f 1,3- |
    # remove commas, which aren't always present within a column
    tr -d ',' |
      awk '
      BEGIN { FS = OFS = "\t" }
  
      FILENAME == ARGV[ARGC-1] {
          if ($1 != orthogroup) 
          if (orthogroup)
              close (orthogroup "_nuc.fa")
              orthogroup = $1
  
          for (col = 2; col <= NF; col++) {
              split($(col), genes, / /)
              for (gene in genes)
                  if (genes[gene] && genes[gene] != "NA") {
                      species_gene = species_col[col-1] "_" genes[gene]
                      print ">" species_gene > orthogroup "_nuc.fa"
                      print fasta[species_gene] > orthogroup "_nuc.fa"
                  }
          }
          next
      }
  
      /^>/ {
         if (! count || species_col[count] != species)
             species_col[++count] = species
             
         if (match($0, /protein_id=[^\]]+/)) # NCBI
             seqid = species "_" substr($0, RSTART+length("protein_id="), RLENGTH-length("protein_id="))
         else # comparative Augustus
             seqid = species "_" substr($1, 2, index($1, " ")-2)
         skip = (seqid in fasta) ? 1 : 0
         if (skip)
             print "WARNING: skipping duplicate gene id: " seqid | "cat 1>&2"
         next
      }
  
      !skip { fasta[seqid] = fasta[seqid] $0 }' \
      species=anaPla /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/anaPla/anaPla_cds_from_genomic.fna \
      species=ansBra /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/ansBra.cds.fa \
      species=ansCyg /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/ansCyg/ansCyg_cds_from_genomic.fna \
      species=ansInd /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/ansInd.cds.fa \
      species=braCan /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/braCan.cds.fa \
      species=colVir /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/colVir/colVir_cds_from_genomic.fna \
      species=cotJap /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/cotJap/cotJap_cds_from_genomic.fna \
      species=galGal /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/galGal/galGal_cds_from_genomic.fna \
      species=hetAtr /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/hetAtr.cds.fa \
      species=netAur /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/netAur.cds.fa \
      species=numMel /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/ncbi_data/numMel/numMel_cds_from_genomic.fna \
      species=oxyJam /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/oxyJam.cds.fa \
      species=stiNae /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/stiNae.cds.fa \
      species=syrMik /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/syrMik.cds.fa \
      species=tymCupPin /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/compAug_data/tymCupPin.cds.fa \
      -
