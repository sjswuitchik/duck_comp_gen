#!/usr/bin/awk -f

# usage: gzip -dc infile.bg.gz | ./sum_cov.awk 

BEGIN { 
    FS = OFS = "\t"
    print "chrom", "start", "end", "weighted.cov" > "stiNae_coverage_sites_low.bed"
    print "chrom", "start", "end", "weighted.cov" > "stiNae_coverage_sites_high.bed"
    print "chrom", "start", "end", "weighted.cov" > "stiNae_coverage_sites_clean.bed"
}
{
 
    int_len = $3 - $2
    mean = 33.87018351181389
    weighted_cov = $4 / int_len
    if (weighted_cov < 0.5*mean)
        print $1, $2, $3, weighted_cov > "stiNae_coverage_sites_low.bed"
    else if (weighted_cov > 2.0*mean)
        print $1, $2, $3, weighted_cov > "stiNae_coverage_sites_high.bed"
    else
        print $1, $2, $3, weighted_cov > "stiNae_coverage_sites_clean.bed"
}
