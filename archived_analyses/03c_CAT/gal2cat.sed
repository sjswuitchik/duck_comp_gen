#!/usr/bin/sed -f
# USAGE
#     gal2cat.sed input.gff > output.cat.gff

# for each "mRNA" and "transcript" feature
/	[mt][Rr][Na][An]/ {
    # convert "gene" attribute to "gene_name" & "gene_id"
    s/gene=\([^;]*\)/gene_name=\1;gene_id=\1/
    # copy "transcript_id" attribute value to new "transcript_name" attribute
    s/transcript_id=\([^;]*\)/transcript_id=\1;transcript_name=\1/
    # add gene_biotype & transcript_biotype attributes
    s/$/;gene_biotype=protein_coding;transcript_biotype=protein_coding/
}
