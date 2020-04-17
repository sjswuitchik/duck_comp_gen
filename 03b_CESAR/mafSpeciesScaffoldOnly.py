#!/usr/bin/env python3
"""A short script to remove version details from src columns in s rows in .maf files."""

# import os
# import sys
import fileinput

def sline_drop_version(line, sep="\t"):
    """Drop version info from src column in s lines."""
    cols_lst = line.split()
    spec_scaf_lst = cols_lst[1].split(".")[:2]
    spec_scaf = ".".join(spec_scaf_lst)

    replaceline = sep.join(["s", spec_scaf] + cols_lst[2:])
    return replaceline


with fileinput.input() as infiles:
    for line in infiles:
        if line.startswith("s"):
            replace = sline_drop_version(line)
            print(replace)
        else:
            print(line, end='')
