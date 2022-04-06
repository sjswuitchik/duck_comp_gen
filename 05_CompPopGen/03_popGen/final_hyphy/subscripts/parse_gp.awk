#!/usr/bin/awk -f
$1 == "LOCUS" { printf("%s", $2) }
$1 ~ /^\/(gene=|db_xref="GeneID)/ { split($1, F, "\""); printf("\t%s", F[2]) }
$1 == "ORIGIN" { printf("\n") }
