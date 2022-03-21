#!/usr/bin/awk -f
# USAGE
#     species_pval.awk absrel_output_clean.csv

BEGIN { 
    FS = OFS = ","
}

NR == 1 {
    print $1,$2,$3,$4,$5,"species:pval"
    next
}

!$6 { print; next }

$6 {
    nf = split($6, species, ";")
    split($7, pvals, ";")
    printf("%s,%s,%s,%s,%s,", $1,$2,$3,$4,$5)
    for (i = 1; i <= nf; i++)
        printf("%s:%s;", species[i], pvals[i])
    printf("\n")
}
