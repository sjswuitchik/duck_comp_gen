#!/bin/perl

use File::Basename;
use strict;

my $in_fasta = $ARGV[0];
my %fa = ();
my $defline;
my $seq;

open IN, "<$in_fasta";
while(my $line=<IN>) {
    chomp $line;
    if($line=~/^>(.*)$/) {
       $defline = $1;
    } 
    else {
	$seq = uc($line);
 	if(!defined($fa{$defline})) {
	    $fa{$defline}{SEQ_LIST}[0]=$seq;
            $fa{$defline}{DUP_COUNT}=1;
        }
        else {
            push @{$fa{$defline}{SEQ_LIST}}, $seq;
	    $fa{$defline}{DUP_COUNT}++;
        }
    }    
}
close IN;

# print out result:
my @suffix_list = (".fa");

my $fa_new = basename($in_fasta, @suffix_list) . "_uniq.fa";
my $fa_dup = basename($in_fasta, @suffix_list) . "_dup.fa";
my $seq_cnt = 0;

open OUT1, ">$fa_new";
open OUT2, ">$fa_dup";

for my $acc (keys %fa) {
    if ($fa{$acc}{DUP_COUNT}==1) {
       print OUT1  ">" . $acc . "\n"; 
       print OUT1 $fa{$acc}{SEQ_LIST}[0] . "\n";
       $seq_cnt++;
    }
    else {
	for my $i (1..$fa{$acc}{DUP_COUNT}) {    
           print OUT2  ">DUP=$fa{$acc}{DUP_COUNT}|SEQ_ID=$i|" . $acc . "\n"; 
	   print OUT2 $fa{$acc}{SEQ_LIST}[$i-1] . "\n";
        } # end for
       $seq_cnt += $fa{$acc}{DUP_COUNT};
    }
}

close OUT1;
close OUT2;

print "total sequences processed = $seq_cnt\n";


