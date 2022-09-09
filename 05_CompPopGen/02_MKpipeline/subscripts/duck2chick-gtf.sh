#!/bin/bash
set -o errexit -o pipefail
# translate into two-column format (DUCK_TRANSCRIPT_ID CHICKEN_TRANSCRIPT_ID)
awk '
BEGIN { FS = OFS = "\t" }
{
    gsub(/,/, "")
    split($1, duck, / /)
    split($2, chicken, / /)
    for (d in duck)
        for (c in chicken)
            print duck[d] "\t" chicken[c]
}' <(tr -d '\r' < hetGal_trans.tsv) |
  # sort by chicken transcript ID
  sort -k 2,2 |
    # join by chicken transcript ID,
    # output 2-column DUCK_TRANSCRIPT_ID CHICKEN_GENE_ID
    join -t $'\t' -1 2 -o 1.1,2.2 - <(sort -k 1,1 transGene_final.tsv) |
     # sort by DUCK_TRANSCRIPT_ID, CHICKEN_GENE_ID, and uniqify output
     sort -u -k 1,1 -k 2,2 |
       # print 2-column: DUCK_TRANSCRIPT_ID CHICKEN_GENE_ID_LIST
       awk '$1 != prev_duck_transcript_id { 
                if (prev_duck_transcript_id)
                    printf("\n")
                printf("%s\t%s", $1, $2)
            }
            $1 == prev_duck_transcript_id { printf(",%s", $2) }
            { prev_duck_transcript_id = $1 }
            END { printf("\n") }' |
         # Augment original GTF with 10th column containing gene list
         awk 'BEGIN { FS = OFS = "\t" }
              NR == FNR { chicken_gene_list[$1] = $2 }
              NR != FNR {
                  if ($3 == "transcript") {
                      transcript_id = $9
                      orig_gene_id = substr(transcript_id, 1, index(transcript_id, ".")-1)
                      gene_id = chicken_gene_list[$9]
                      gsub(/,/, "-", gene_id) # concatenate multiple gene IDs w/ hyphens
                      $9 = gene_id ".t1"
                  } else if (substr($1,1,1) != "#")
                      gsub(orig_gene_id, gene_id, $9)
                  print
              }' - hetAtr.gff
