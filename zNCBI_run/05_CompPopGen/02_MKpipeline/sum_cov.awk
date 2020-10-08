#!/usr/bin/awk -f

# usage: gzip -dc infile.bg.gz | ./sum_cov.awk 

BEGIN { 
    FS = OFS = "\t"
    print "chrom", "start", "end", "weighted.cov" | "gzip > coverage_sites_low.bed.gz"
    print "chrom", "start", "end", "weighted.cov" | "gzip > coverage_sites_high.bed.gz"
    print "chrom", "start", "end", "weighted.cov" | "gzip > coverage_sites_clean.bed.gz"
}
{
    tot_cov = 0 
    for (ind = 4; ind <= 14; ind++)
        tot_cov += $ind
    int_len = $3 - $2
    weighted_cov = tot_cov / int_len
    if (weighted_cov < 0.5)
        print $1, $2, $3, weighted_cov | "gzip > coverage_sites_low.bed.gz"
    else if (weighted_cov > 2.0)
        print $1, $2, $3, weighted_cov | "gzip > coverage_sites_high.bed.gz"
    else
        print $1, $2, $3, weighted_cov | "gzip > coverage_sites_clean.bed.gz"
}
