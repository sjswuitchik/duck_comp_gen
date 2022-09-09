#!/usr/bin/env python3
############################################################
# Reads results from Hyphy json output
############################################################

import sys, os, argparse, lib.hpcore as hpcore

############################################################
# Options

parser = argparse.ArgumentParser(description="Parse Hyphy json output");
parser.add_argument("-i", dest="input", help="Directory containing hyphy json output files.", default=False);
parser.add_argument("-m", dest="model", help="The Hyphy model that was used to generate the files in -i. Default: mg94, mg94-local, anc-recon, busted, fel, fubar, absrel, slac", default="mg94");
parser.add_argument("-d", dest="meta", help="A file containing metadata. Tab delimited with columns: id, feature type, chromosome, start coord, end coord, strand. Directories in -i must be formatted <id>-*", default=False);
parser.add_argument("-o", dest="output", help="An output .csv file.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
args = parser.parse_args();

if not args.input or not os.path.isdir(args.input):
    sys.exit(" * Error 1: Please provide a valid input directory (-i).");

if args.model not in ["mg94", "mg94-local", "anc-recon", "busted", "fel", "fubar", "absrel", "slac"]:
    sys.exit(" * Error 2: Model (-m) must be one of: mg94, mg94-local, anc-recon, busted, fel, fubar, absrel, slac");

if args.meta and not os.path.isfile(args.meta):
    sys.exit(" * Error 3: Cannot find meta data file: " + args.meta);

if not args.output:
    sys.exit(" * Error 2: Please provide the name of an output file (-o).")

if os.path.isfile(args.output) and not args.overwrite:
    sys.exit( " * Error 3: Output file (-o) already exists! Explicity specify --overwrite to overwrite it.");

site_models =  ["fubar"]
if args.model in site_models:
    if args.output.endswith(".csv"):
        sites_file = args.output.replace(".csv", ".sites.csv");
    else:
        sites_file = args.output + ".sites";

pad = 25;

with open(args.output, "w") as outfile:
    hpcore.runTime("# HyPhy output parser", outfile);
    hpcore.PWS("# IO OPTIONS", outfile);
    hpcore.PWS(hpcore.spacedOut("# Input directory:", pad) + args.input, outfile);
    hpcore.PWS(hpcore.spacedOut("# Hyphy model:", pad) + args.model, outfile);
    if args.meta:
        hpcore.PWS(hpcore.spacedOut("# Metadata file:", pad) + args.meta, outfile);
    hpcore.PWS(hpcore.spacedOut("# Output file:", pad) + args.output, outfile);
    if args.model in site_models:
         hpcore.PWS(hpcore.spacedOut("# Sites file:", pad) + sites_file, outfile);
    if args.overwrite:
        hpcore.PWS(hpcore.spacedOut("# --overwrite set:", pad) + "Overwriting previous output file.", outfile);
    hpcore.PWS("# ----------------", outfile);

    features = False;
    if args.meta:
        hpcore.PWS("# " + hpcore.getDateTime() + " Reading metadata file: " + args.meta, outfile);
        features = hpcore.readMeta(args.meta);
        hpcore.PWS(hpcore.spacedOut("# Features read:        ", pad) + str(len(features)), outfile);               
        hpcore.PWS("# ----------------", outfile);
    # Read the feature metadata.

    hpcore.PWS("# " + hpcore.getDateTime() + " Begin parsing Hyphy output...", outfile);
    if args.model == "mg94":
        import lib.mg94 as mg94;
        mg94.parse(args.input, features, outfile, pad);

    if args.model == "mg94-local":
        import lib.mg94local as mg94local;
        mg94local.parse(args.input, features, outfile, pad);

    if args.model == "anc-recon":
        import lib.ancrecon as ancrecon;
        ancrecon.parse(args.input, features, outfile, pad);

    if args.model == "busted":
        import lib.busted as busted;
        busted.parse(args.input, features, outfile, pad);

    if args.model == "fel":
        import lib.fel as fel;
        fel.parse(args.input, features, outfile, pad);

    if args.model == "fubar":
        import lib.fubar as fubar;
        fubar.parse(args.input, features, outfile, sites_file, pad);

    if args.model == "absrel":
        import lib.absrel as absrel;
        absrel.parse(args.input, features, outfile, pad);

    if args.model == "slac":
        import lib.slac as slac;
        slac.parse(args.input, features, outfile, pad);

    if args.model == "relax":
        import lib.relax as relax;
        relax.parse(args.input, features, outfile, pad);
    # Load the library for the model used and pass everything to it.
