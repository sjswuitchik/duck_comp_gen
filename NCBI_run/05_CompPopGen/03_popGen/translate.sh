#!/bin/bash
set -o errexit -o nounset
for fasta in single_orthos/*.fa
do
  awk 'NR == FNR {species[$2] = $1}
       NR != FNR {
           if (/^>/) {
               sp = species[substr($1,2)]
               if (sp == "") {
                   print FILENAME ":" NR " - Species not found for seqid: " substr($1, 2) | "cat 1>&2"
                   exit(1)
               }
               print ">" sp
           } else
               print
       }' <(singularity exec /cvmfs/singularity.galaxyproject.org/n/e/newick_utils:1.6--h516909a_3 nw_labels gene_trees/$(basename ${fasta} .fa)_tree.txt | sed $'s/_translated_/\t/g' ) ${fasta} > translated/$(basename ${fasta})
done
