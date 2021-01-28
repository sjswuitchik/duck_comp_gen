#!/usr/bin/awk -f

# usage: gzip -dc infile.bg.gz | ./sum_cov.awk 

BEGIN { 
    FS = OFS = "\t"
    print "chrom", "start", "end", "weighted.cov" > "hetAtr_coverage_sites_low.bed"
    print "chrom", "start", "end", "weighted.cov" > "hetAtr_coverage_sites_high.bed"
    print "chrom", "start", "end", "weighted.cov" > "hetAtr_coverage_sites_clean.bed"
}
{
    tot_cov = 0 
    for (ind = 4; ind <= NF; ind++)
        tot_cov += $ind
    int_len = $3 - $2
    weighted_cov = tot_cov / int_len
    if (weighted_cov < 0.5)
        print $1, $2, $3, weighted_cov > "hetAtr_coverage_sites_low.bed"
    else if (weighted_cov > 2.0)
        print $1, $2, $3, weighted_cov > "hetAtr_coverage_sites_high.bed"
    else
        print $1, $2, $3, weighted_cov > "hetAtr_coverage_sites_clean.bed"
}
