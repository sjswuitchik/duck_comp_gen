#!/usr/bin/python
############################################################
# For spider genomes, 01.2022
# Trims and filters codon alignments
## from https://github.com/gwct/murine-discordance/blob/main/scripts/03-selection-tests/08_aln_filter.py
############################################################

import sys, os, argparse
sys.path.append("/home/gt156213e/bin/core/python/lib/");
import core as CORE
import seqparse as SEQ
#from Bio.Seq import Seq

############################################################
# Functions

def ntToCodon(nt_seq):
# This function takes a string of sequence and parses it into a list of codons.
# String is required to be divisible by 3.

    assert len(nt_seq) % 3 == 0, "\nOUT OF FRAME NUCLEOTIDE SEQUENCE! " + str(len(nt_seq));
    # Check that sequence is in frame.

    codon_seq = [(nt_seq[i:i+3]) for i in range(0, len(nt_seq), 3)];
    # Get chunks of 3 characters into a list.

    return codon_seq;

#########################

def countNongapLength(codon_seqs, seq_filter):
# This function goes through every sequence in an alignment and calculates
# the average length of each sequence excluding gaps.

    len_sum, num_seqs, num_gappy_seqs, gappy_seqs = 0, 0, 0, [];
    for seq in codon_seqs:
        num_seqs += 1;
        full_len = len(codon_seqs[seq]);
        nogap_len = len( [ codon for codon in codon_seqs[seq] if codon != "---" ] );
        len_sum += nogap_len;

        if 1 - (nogap_len / full_len) > seq_filter:
            num_gappy_seqs += 1;
            gappy_seqs.append(seq);
        # If the number of gaps in the sequence (calculated as 1 - the fraction of non-gap length to full length)
        # is above some threshold, mark it as gappy

    return len_sum / num_seqs, num_gappy_seqs, gappy_seqs;

#########################

def countUniqIdentSeqs(codon_seqs, gappy_seqs, post_stop_seqs):
# This function goes through every sequence in an alignment and counts how 
# many sequences are unique or identical.

    uniq_seqs, ident_seqs, found = 0, 0, [];
    #codon_seq_list = list(codon_seqs.values());
    codon_seq_list = [ codon_seqs[seq] for seq in codon_seqs if seq not in gappy_seqs + post_stop_seqs ]

    for seq in codon_seq_list:
        if codon_seq_list.count(seq) == 1:
            uniq_seqs += 1;
        if codon_seq_list.count(seq) != 1 and seq not in found:
            ident_seqs += 1;
            found.append(seq);

    return uniq_seqs, ident_seqs;

#########################

def siteCount(nt_seqs, codon_seqs, aln_len_codons, aln_name):
# This function goes through every site in an alignment to check for invariant sites, gappy sites, 
# and stop codons.

    #stop_codons = ["TAG", "TAA", "TAR", "TGA", "TRA"];
    stop_codons = ["TAG", "TAA", "TGA", "UAG", "UAA", "UGA"];
    # Possible stop codons with IUPAC ambiguities

    invar_sites, gap_sites, stop_codon_count, high_gap_sites = 0, 0, 0, 0;
    informative_sites_codon = 0;
    # Counts

    codon_seq_list = list(codon_seqs.values());
    # if "ENSMUST00000040971" in aln_name:
    #     print(codon_seq_list);
    #     sys.exit();
    for i in range(aln_len_codons):
    # Loop over every site

        site = [];
        for j in range(len(codon_seq_list)):
            site.append(codon_seq_list[j][i]);
        # Get the current (ith) site from every sequence (j)

        if site.count(site[0]) == len(site):
            invar_sites += 1;
        # If all the codons at the site match the first one, it is invariant

        num_gaps = site.count("---");
        # Count the number of gaps at the current site

        allele_counts = { allele : site.count(allele) for allele in site if not any(missing_char in allele for missing_char in ["-", "N", "X"]) };
        # Count the occurrence of each allele in the site

        if len(allele_counts) > 1:
            multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
            # Count the number of allele present in at least 2 species

            if len(multi_allele_counts) >= 2:
                informative_sites_codon += 1;
            # If 2 or more alleles are present in 2 or more species, this site is informative

        if num_gaps > 1:
            gap_sites += 1;
            # Increment by one if there is at least one gap

            if (num_gaps / len(site)) > 0.2:
                high_gap_sites += 1;
            # Count if the number of gaps at this site is above some threshold
        
        if i != aln_len_codons-1:
            stop_codon_count += len( [ codon for codon in site if codon in stop_codons ] );
        # Count the number of stop codons at the site
    ## End site loop

    stop_codon_seqs = [];
    for sample in codon_seqs:
        for codon in codon_seqs[sample][:-1]:
            if codon in stop_codons:
                #print(sample, codon);
                stop_codon_seqs.append(sample);
                break;
    # Get a list of which sequences have stop codons
    ## TODO: add exception for last codon?

    aln_len_nt = aln_len_codons * 3;
    informative_sites_nt = 0;
    nt_seq_list = list(nt_seqs.values());
    for i in range(aln_len_nt):

        site = [];
        for j in range(len(nt_seq_list)):
            site.append(nt_seq_list[j][i]);

        allele_counts = { allele : site.count(allele) for allele in site if allele not in ["-", "N", "X"] };
        # Count the occurrence of each allele in the site

        if len(allele_counts) > 1:
            multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
            # Count the number of allele present in at least 2 species

            if len(multi_allele_counts) >= 2:
                informative_sites_nt += 1;
            # If 2 or more alleles are present in 2 or more species, this site is informative
    ## TODO: Count parsimony informative sites

    return invar_sites, gap_sites, stop_codon_count, high_gap_sites, stop_codon_seqs, informative_sites_nt, informative_sites_codon;

#########################

def codonWindowFilter(codon_seqs, num_samples, cds_len, site_filter):
# This function takes a sliding window along a codon alignment and filters windows
# where 2 or more codons have 2 or more gaps.

    exclude_codons = [];
    # The indices of the sites to filter

    i = 0;
    while i < cds_len-3:
    # Loop over every codon in the alignment, stopping at the last possible 3 codon window

        seq_exclude = 0;
        # The number of sequences at the current codon window that are too gappy (at least 2 gaps in at least 2 codons)

        for seq in codon_seqs:
        # Go over every sequence in the alignment

            c1 = codon_seqs[seq][i];
            c2 = codon_seqs[seq][i+1];
            c3 = codon_seqs[seq][i+2]; 
            # For every sequence get the current codon and the next two codons

            gappy_c = [ c for c in [c1, c2, c3] if c.count("-") >= 2 ];
            # Of the 3 codons in the current window, get the ones with 2 or more gaps

            if len(gappy_c) >= 2:
                seq_exclude += 1;
            # If there are 2 or more gappy codons in the window in the current sequence, add it to the
            # list of exclude sequences
        ## End seq loop

        if seq_exclude / num_samples > site_filter:
            exclude_codons += [i, i+1, i+2];
        # If the total number of sequences excluded for being gappy at these sites is over some threshold,
        # add all 3 codons to the list of excluded sites

        i += 1;
    ## End codon loop

    for seq in codon_seqs:
        codon_seqs[seq] = [ codon_seqs[seq][i] for i in range(cds_len) if i not in exclude_codons ];
    # From every sequence, remove the codons determined to be gappy in too many sequences above

    return codon_seqs, len(list(set(exclude_codons)));

#########################

def writeAln(seqs, out_file):
# This function writes sequences to a file in FASTA format.

    with open(out_file, "w") as outfile:
        for seq in seqs:
            outfile.write(seq + "\n");
            outfile.write(seqs[seq] + "\n");

############################################################
# Options

parser = argparse.ArgumentParser(description="Codon alignment check filter");
parser.add_argument("-i", dest="input", help="Directory of CDS alignments in fasta format.", default=False);
parser.add_argument("-n", dest="num_spec_filter", help="Alignments with sequences from fewer than this many taxa after filtering will be removed. Default: 4", type=int, default=4);
parser.add_argument("-d", dest="header_delim", help="FASTA headers should be in the form of <gene id>_<species id>. Use this to change the delimiter from _ to another character. For ' ' (space), enter 'space'.", default="_");
parser.add_argument("-s", dest="seq_filter", help="The gappy sequence filter threshold. Sequences that are greater than this percent gappy will be removed. Default: 20", type=int, default=20);
parser.add_argument("-c", dest="site_filter", help="The codon threshold. Sites that have greater than this percent sequences fail the codon window filter will be removed. Default: 50", type=int, default=50);
parser.add_argument("-o", dest="output", help="Desired output directory for filtered CDS alignments and translated AA alignments.", default=False);
# User params

parser.add_argument("--count", dest="count_only", help="Set this option to just provide the log file with counts/stats. Will not filter or write new sequences", action="store_true", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# User options

args = parser.parse_args();

## Option parsing
##########################

if not args.input or not os.path.isdir(args.input):
    sys.exit( " * Error 1: An input directory with aligned CDS sequences must be defined with -i.");
args.input = os.path.abspath(args.input);
# Check input directory path

if not args.num_spec_filter:
    sys.exit( " * Error 2: Please provide the previous presence filter (-f).");
# Check the spec filter

if args.header_delim == "space":
    args.header_delim = " ";
# Parse the FASTA header delimiter for species labels

if args.seq_filter < 0 or args.seq_filter > 100:
    sys.exit( " * Error 3: Sequence gap threshold (-g) must be between 0 and 100.");
# Check the gappy sequence filter

if args.site_filter < 0 or args.site_filter > 100:
    sys.exit( " * Error 4: Codon window gap threshold (-c) must be between 0 and 100.");
# Check the gappy codon window/site filter

if not args.count_only and not args.output:
    sys.exit( " * Error 5: An output directory must be defined with -o.");
# Make sure an output directory is provided

## Input options and params
##########################

args.output = os.path.abspath(args.output);
if args.output[-1] == "/":
    args.output = args.output[:-1];
args.output += "-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter);
# Adjust the output directory name to add the filter thresholds.

if os.path.isdir(args.output) and not args.overwrite:
    sys.exit( " * Error 6: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# Check the overwrite option

nt_outdir = os.path.join(args.output, "cds");
aa_outdir = os.path.join(args.output, "pep");
if not os.path.isdir(args.output):
    os.system("mkdir " + args.output);
if not args.count_only:
    os.system("mkdir " + nt_outdir);
    os.system("mkdir " + aa_outdir);
# Create the output directory and sub-directories for both alignment types

spec_file = os.path.join(args.output, "spec-stats-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# Main output file with stats per sequence

log_file = os.path.join(args.output, "aln-stats-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".log");
# Log file

passed_loci_file = os.path.join(args.output, "aln-filter-passed-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".txt");

rm_too_few_file = os.path.join(args.output, "too-few-species-filtered-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# A file that contains a list of alignments that are removed because they have too few species remaining after filtering

rm_too_short_file = os.path.join(args.output, "too-short-filtered-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# A file that contains a list of alignments that are removed because they are too short filtering

rm_stop_file = os.path.join(args.output, "stop-codon-filtered-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# A file that contains a list of sequences that are removed because they have premature stop codons

rm_gappy_file = os.path.join(args.output, "gappy-seqs-filtered-spec" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# A file that contains a list of sequences that are removed because they exceed the gappy sequence threshold

#rm_protein_file = os.path.join(args.output, "gappy-proteins-f" + str(args.num_spec_filter) + "-seq" + str(args.seq_filter) + "-site" + str(args.site_filter) + ".tab");
# A file that contains a list of sequences that are removed because they have premature stop codons

# Output options and files
##########################

pad = 26
# Job vars

with open(log_file, "w") as logfile, open(spec_file, "w") as specfile:
# Open main output and log files

    CORE.runTime("# CDS alignment filter", logfile);
    CORE.PWS("# IO OPTIONS", logfile);
    CORE.PWS(CORE.spacedOut("# Input CDS directory:", pad) + args.input, logfile);
    CORE.PWS(CORE.spacedOut("# Sequence gappiness threshold:", pad) + str(args.seq_filter), logfile);
    CORE.PWS(CORE.spacedOut("# Codon site gappiness threshold:", pad) + str(args.site_filter), logfile);
    CORE.PWS(CORE.spacedOut("# Min sequences required after filters:", pad) + str(args.num_spec_filter), logfile);
    CORE.PWS(CORE.spacedOut("# Output directory:", pad) + args.output, logfile);
    if args.overwrite:
        CORE.PWS(CORE.spacedOut("# --overwrite set:", pad) + "Overwriting previous files in output directory.", logfile);
    if args.count_only:
        CORE.PWS(CORE.spacedOut("# --count set:", pad) + "Will not output sequences.", logfile);
    CORE.PWS(CORE.spacedOut("# Log file:", pad) + log_file, logfile);
    CORE.PWS("# ----------------", logfile);
    # Report options at runtime and in the log
    
    args.seq_filter = args.seq_filter / 100;
    args.site_filter = args.site_filter / 100;
    # Adjust the filters to be proportions

    ## Reporting run-time info for records.
    ########################## 

    CORE.PWS("# " + CORE.getDateTime() + " Beginning filter.", logfile);
    # Status update

    witten_loci, rm_too_few_spec, rm_too_short, rm_stop_codons, rm_gappy, rm_protein_gappy = [], [], [], [], [], [];
    # Lists of aligns/sequences removed by each filter

    pre_single_copy_alns = [];
    post_single_copy_alns = [];
    num_spec = 7;

    pre_alns, pre_seqs, post_alns, post_seqs = 0, 0, 0, 0;
    # Counts of aligns and sequences before and after filtering

    aln_stats = ["num.seqs", "num.species", "codon.aln.length", "avg.nongap.length", "uniq.seqs", "ident.seqs", "gappy.seqs", "invariant.sites", "informative.codon.sites", "informative.nt.sites", "stop.codons", 
                "percent.sites.with.gap", "gappy.sites"];
    # Stats to collect

    aln_headers = ["align"] + [ "pre." + s for s in aln_stats ];
    if not args.count_only:
        aln_headers += ["sites.filtered"] + [ "post." + s for s in aln_stats ];
    # Headers for the output file including some that need to appear only once, plus some the stats that appear both pre and post filter

    CORE.PWS("\t".join(aln_headers), logfile);
    # Write the headers to the log file

    ##########################

    spec_headers = ["spec"] + ["num.alns", "num.gappy", "num.stop"];
    specfile.write("\t".join(spec_headers) + "\n");
    # Write the headers to the species stats output file

    spec_stats = {};
    # The species counts dict

    ##########################
    
    aln_files = [ f for f in os.listdir(args.input) if os.path.isfile(os.path.join(args.input, f)) and f.endswith(".fa") ];
    num_alns = len(aln_files);
    num_alns_str = str(num_alns);
    num_alns_len = len(num_alns_str);
    # Read alignment directory names from input directory

    ##########################

    written, num_short, num_high_ident, num_gappy, aln_prem_stop, num_no_info, num_stoppy = 0.0,0.0,0.0,0.0,0.0,0.0,0.0;
    # Some count variables for all aligns

    first_aln, counter, skipped = True, 0, 0;
    # Loop tracking variables

    for aln in aln_files:
    # Loop over every alignment

        # if "ENSMUST00000024939" not in aln:
        #     continue;

        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_alns_len:
                counter_str = "0" + counter_str;
            print ("> " + CORE.getDateTime() + " " + counter_str + " / " + num_alns_str);
        counter += 1;
        # Loop progress   

        pre_alns += 1;
        # Increment the alignment counter

        cur_out = { h : "NA" for h in aln_headers };
        cur_out["align"] = aln;
        # Initialize the current output dictionary with NAs for most columns

        cur_infile = os.path.join(args.input, aln);
        # The alignment file

        ##########################

        if not args.count_only:
            cur_nt_outfile = os.path.join(nt_outdir, aln.replace(".fa", ".filter.NT.fa"));
            cur_aa_outfile = os.path.join(aa_outdir, aln.replace(".fa", ".filter.AA.fa"));
        # Get the current in and output files

        ##########################

        seqs_orig = SEQ.fastaGetDict(cur_infile);
        seqs = { t : seqs_orig[t].upper().replace("!", "-") for t in seqs_orig };
        # Read the sequences and convert all chars to upper case and remove MACSE's exclamation marks that indicate frameshifts

        cur_seq_names = list(seqs.keys());
        pre_num_seq = len(cur_seq_names);
        pre_seqs += pre_num_seq;
        # Counts the number of sequences in the alignment before filtering

        cur_spec_names = set([ s.split(args.header_delim)[0] for s in cur_seq_names ]);
        pre_num_spec = len(cur_spec_names);
        if pre_num_spec == num_spec and pre_num_seq == num_spec:
            pre_single_copy_alns.append(aln);
        # Counts the number of species represented in the alignment before filtering

        for spec in cur_spec_names:
            if spec not in spec_stats:
                spec_stats[spec] = { col : 0 for col in spec_headers if col != "spec" };
            spec_stats[spec]['num.alns'] += 1;
        # Count the species in the alignment in the species dict and initialize if it is the first time this species is seen

        codons = {};
        for seq in seqs:    
            codons[seq] = ntToCodon(seqs[seq]);
        # Convert from NT strings to codon lists

        ##########################

        cur_out["pre.num.seqs"] = str(pre_num_seq);
        # Count the number of samples in the alignment

        cur_out["pre.num.species"] = str(pre_num_spec);

        cur_out["pre.codon.aln.length"] = len(codons[cur_seq_names[0]]);
        # Count the overall length of the alignment

        cur_out["pre.avg.nongap.length"], cur_out["pre.gappy.seqs"], gappy_seqs = countNongapLength(codons, args.seq_filter);
        # Calculate the average length of the sequences in the alignment without gaps and the number of sequences that are over the gap threshold

        cur_out["pre.uniq.seqs"], cur_out["pre.ident.seqs"] = countUniqIdentSeqs(seqs, [], []);
        # Counts the number of unique sequences in an alignment (uniq + ident = total unique sequences)

        cur_out["pre.invariant.sites"], cur_out["pre.percent.sites.with.gap"], cur_out["pre.stop.codons"], cur_out["pre.gappy.sites"], pre_stop_seqs, cur_out["pre.informative.nt.sites"], cur_out["pre.informative.codon.sites"] = siteCount(seqs, codons, cur_out["pre.codon.aln.length"], aln);
        # Count invariant, informative, gaps, and stop codons at each site in the alignment

        ##########################

        if not args.count_only:
        # Do some filtering if we aren't just counting the input sequences

            f_codons, cur_out["sites.filtered"] = codonWindowFilter(codons, pre_num_seq, cur_out["pre.codon.aln.length"], args.site_filter);
            # Remove sites that are gappy in many sequences (2 gaps in at least 2 codons in overlapping 3 codon windows in at least 50% of sequences)

            f_seqs = { seq_id : "".join(f_codons[seq_id]) for seq_id in f_codons };
            # Convert from codon lists to sequence strings

            cur_out["post.avg.nongap.length"], cur_out["post.gappy.seqs"], gappy_seqs = countNongapLength(f_codons, args.seq_filter);
            # Calculate the average length of the sequences in the alignment without gaps and the number of sequences that are over the gap threshold AFTER FILTERING CODON SITES

            cur_out["post.codon.aln.length"] = len(f_codons[cur_seq_names[0]]);
             # Count the overall length of the alignment AFTER FILTERING CODON SITES

            short_aln = False;
            if cur_out["post.codon.aln.length"] < 33:
                rm_too_short.append(aln);
                short_aln = True;
                num_short += 1;

            cur_out["post.invariant.sites"], cur_out["post.percent.sites.with.gap"], cur_out["post.stop.codons"], cur_out["post.gappy.sites"], post_stop_seqs, cur_out["post.informative.nt.sites"], cur_out["post.informative.codon.sites"] = siteCount(f_seqs, f_codons, cur_out["post.codon.aln.length"], aln);
            # Count invariant, informative, gaps, and stop codons at each site in the alignment AFTER FILTERING CODON SITES

            cur_out["post.uniq.seqs"], cur_out["post.ident.seqs"] = countUniqIdentSeqs(f_codons, gappy_seqs, post_stop_seqs);
            # Counts the number of unique sequences in an alignment (uniq + ident = total unique sequences) AFTER FILTERING CODON SITES

            ##########################

            for seq in gappy_seqs:
                rm_gappy.append(aln + "\t" + seq);
                gappy_spec = seq.split(args.header_delim)[0];
                spec_stats[gappy_spec]["num.gappy"] += 1;
            # Add the seqs removed for being too gappy to the gappy seq list

            for seq in post_stop_seqs:
                rm_stop_codons.append(aln + "\t" + seq);
                stop_spec = seq.split(args.header_delim)[0];
                spec_stats[stop_spec]["num.gappy"] += 1;
            # Add the seqs removed for being having stop codons to the stop codon seq list

            ##########################

            stop_codons = ["TAG", "TAA", "TGA", "UAG", "UAA", "UGA"];
            #print(f_codons);
            if any(f_codons[seq][-1] in stop_codons for seq in f_codons):
                f_codons = { seq : f_codons[seq][:-1] for seq in f_codons };
                f_seqs = { seq : f_seqs[seq][:-3] for seq in f_seqs };
            #print(f_codons);
            # Remove the trailing stop codons...

            final_codons = { seq : f_codons[seq] for seq in f_codons if seq not in gappy_seqs and seq not in post_stop_seqs };
            final_seqs = { seq : f_seqs[seq] for seq in final_codons if seq not in gappy_seqs and seq not in post_stop_seqs };
            # Remove sequences from the alignment that are too gappy or have stop codons after codon filtering

            post_num_seqs = len(final_seqs);
            cur_out["post.num.seqs"] = post_num_seqs;
            post_spec_names = set([ s.split(args.header_delim)[0] for s in final_seqs ]);
            post_num_spec = len(post_spec_names);
            if post_num_spec == num_spec and post_num_seqs == num_spec:
                post_single_copy_alns.append(aln);
            cur_out["post.num.species"] = post_num_spec;
            # Count the number of sequences and species in the alignment after filtering

            write_aln = True;
            if cur_out["post.num.species"] < args.num_spec_filter:
                rm_too_few_spec.append(aln);
                write_aln = False;
            elif short_aln:
                write_aln = False;
            else:
                post_alns += 1;
                post_seqs += len(final_seqs);
            # If after filtering there are fewer species left than required by the threshold, add the whole orthogroup to the list of groups with too few species
            # Otherwise increment the counters

            ##########################

            if write_aln:
                writeAln(final_seqs, cur_nt_outfile);
                # NT

                #aas = { sample : str(Seq(f_seqs[sample]).translate()) for sample in f_seqs };
                final_aas = { seq : SEQ.yabt(final_codons[seq]) for seq in final_codons };
                writeAln(final_aas, cur_aa_outfile);
                # AA

                witten_loci.append(aln);
            # Write both the nt and aa sequences if the alignment isn't removed
        ## End filter/write block

        outline = [ str(cur_out[col]) for col in aln_headers ];
        logfile.write("\t".join(outline) + "\n");
        # Write the output per alignment to the log file
    ## End alignment loop

    ##########################

    for spec in spec_stats:
        outline = [spec] + [ str(spec_stats[spec][col]) for col in spec_headers if col != "spec" ];
        specfile.write("\t".join(outline) + "\n");

    ##########################

    CORE.PWS("# ----------------", logfile);
    CORE.PWS("# Pre-filter alignments: " + str(pre_alns), logfile);
    CORE.PWS("# Pre-filter sequences : " + str(pre_seqs), logfile);
    CORE.PWS("# Pre-filter single-copy alignments: " + str(len(pre_single_copy_alns)), logfile);
    # Report the pre-filter counts

    ##########################

    with open(rm_too_few_file, "w") as fewfile:
        for aln in rm_too_few_spec:
            fewfile.write(aln + "\n");
    CORE.PWS("# Alignments removed with too few species: " + str(len(rm_too_few_spec)), logfile);
    # Write the alignments removed with too few species

    ##########################

    with open(rm_too_short_file, "w") as shortfile:
        for aln in rm_too_short:
            shortfile.write(aln + "\n");
    CORE.PWS("# Alignments removed for being too short: " + str(len(rm_too_short)), logfile);
    # Write the alignments removed that are too short

    ##########################

    with open(rm_gappy_file, "w") as gappyfile:
        for seq in rm_gappy:
            gappyfile.write(seq + "\n");
    CORE.PWS("# Samples removed above gappy threshold:    " + str(len(rm_gappy)), logfile);
    # Write the sequence names removed for being too gappy

    ##########################

    with open(rm_stop_file, "w") as stopfile:
        for seq in rm_stop_codons:
            stopfile.write(seq + "\n");
    CORE.PWS("# Samples removed with premature stops:    " + str(len(rm_stop_codons)), logfile);
    # Write the the sequence names removed with premature stop codons

    ##########################

    CORE.PWS("# Post-filter alignments: " + str(post_alns), logfile);
    CORE.PWS("# Post-filter sequences : " + str(post_seqs), logfile);
    CORE.PWS("# Post-filter single-copy alignments: " + str(len(post_single_copy_alns)), logfile);

    ##########################

    with open(passed_loci_file, "w") as passedfile:
        for aln in witten_loci:
            passedfile.write(aln + "\n"); 
    # Write the loci that passed all filters to a file

    ##########################

## Close output and log files

############################################################
