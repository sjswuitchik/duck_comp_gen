#!/usr/bin/awk -f

# usage: gzip -dc hetAtr.merge.bg.gz | awk -v avg=296 -v spp=hetAtr -f sum_cov.awk

BEGIN { 
    FS = OFS = "\t"
    print "chrom", "start", "end", "weighted.cov" > spp"_coverage_sites_low.bed"
    print "chrom", "start", "end", "weighted.cov" > spp"_coverage_sites_high.bed"
    print "chrom", "start", "end", "weighted.cov" > spp"_coverage_sites_clean.bed"
}
{
    mean = avg
    int_len = $3 - $2
    weighted_cov = $4 / int_len
    if (weighted_cov < 0.5*mean)
        print $1, $2, $3, weighted_cov > spp"_coverage_sites_low.bed"
    else if (weighted_cov > 2.0*mean)
        print $1, $2, $3, weighted_cov > spp"_coverage_sites_high.bed"
    else
        print $1, $2, $3, weighted_cov > spp"_coverage_sites_clean.bed"
}
