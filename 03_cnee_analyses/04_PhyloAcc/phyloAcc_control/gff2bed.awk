#
# converts a GFF file to BED format
#
# usage:
#
# $ awk -f gff2bed.awk input.gff > output.bed
#

BEGIN {
    # splitting character is the tab
    FS = OFS = "\t"
}

# skip  headers
/^#/ {
    next;
}

# transform the output
{
    print $1, $4-1, $5, $3, $9 
}
