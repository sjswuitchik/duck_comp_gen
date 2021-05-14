#!/bin/bash
# SYNPSIS
#     augustus-gtf-transcript-to-gene.sh GTF_FILE TRANSCRIPT_GENE_FILE SPECIES_TRANSLATION_FILE 
#
# INPUT FILES
#     GTF_FILE
#         A comparative Augustus GTF file
#
#     TRANSCRIPT_GENE_FILE
#         Two tab-separated columns:
#         1. (reference species) transcript ID
#         2. gene ID
#
#     SPECIES_TRANSLATION_FILE
#         Two tab-separated columns
#         1. Comma-Deparated list of (reference) species transcript ID(s)
#         2. Comma-separated list of Comparative Augustus transcript_id's
#
# STDOUT
#     Input GTF_FILE with 
#
# EXAMPLES
#    $ awk -F '\t' '$3 == "transcript"' ansBra.gff | head -n 1
#    NXHY01000637.1  AUGUSTUS    transcript  35439   36548   .   -   .   jg1.t1
#    $ grep jg1.t1 ansBra_trans.tsv
#    rna-XM_004945648.3, rna-XM_015295016.2   jg1.t1
#    $ grep rna-XM_015295016.2 transGene_final.tsv
#    rna-XM_015295016.2      LOC101751902   
#    $ augustus-gtf-transcript-to-gene.sh ansBra.gff transGene_final.tsv ansBra_trans.tsv
#    ...
#    NXHY01000637.1  AUGUSTUS    transcript  35439   36548   .   -   .   LOC101751902.t1
#    ...
#
# AUTHOR
#     Nathan Weeks <nweeks@g.harvard.edu>

set -o errexit -o pipefail

# translate into two-column format (REF_TRANSCRIPT_ID AUGUSTUS_TRANSCRIPT_ID)
awk '
BEGIN { FS = OFS = "\t" }
{
    gsub(/,/, "")
    split($1, ref_transcript_ids, / /)
    split($2, augustus_transcript_ids, / /)
    for (t in ref_transcript_ids)
        for (a in augustus_transcript_ids)
            print ref_transcript_ids[t] "\t" augustus_transcript_ids[a]
}' <(tr -d '\r' < "${3}") |
  # sort by reference transcript ID
  sort -k 1,1 |
    # join by reference transcript ID,
    # output 2-column AUGUSTUS_TRANSCRIPT_ID, GENE_ID
    join -t $'\t' -o 1.2,2.2 - <(sort -k 1,1 "${2}") |
     # sort by AUGUSTUS_TRANSCRIPT_ID, GENE_ID, and uniqify output
     sort -u -k 1,1 -k 2,2 |
       # print 2-column: AUGUSTUS_TRANSCRIPT_ID GENE_ID_LIST
       awk '$1 != prev_augustus_transcript_id { 
                if (prev_augustus_transcript_id)
                    printf("\n")
                printf("%s\t%s", $1, $2)
            }
            $1 == prev_augustus_transcript_id { printf(",%s", $2) }
            { prev_augustus_transcript_id = $1 }
            END { printf("\n") }' |
         # Augment original GTF with 10th column containing gene list
         awk 'BEGIN { FS = OFS = "\t" }
              NR == FNR { gene_list[$1] = $2 }
              NR != FNR {
                  if ($3 == "transcript") {
                      transcript_id = $9
                      orig_gene_id = substr(transcript_id, 1, index(transcript_id, ".")-1)
                      gene_id = gene_list[$9]
                      if (gene_id) {
                          gsub(/,/, "-", gene_id) # concatenate multiple gene IDs w/ hyphens
                          $9 = gene_id ".t1"
                      }
                  } else if (gene_id && substr($1,1,1) != "#")
                      gsub(orig_gene_id, gene_id, $9)
                  print
              }' - "${1}"
