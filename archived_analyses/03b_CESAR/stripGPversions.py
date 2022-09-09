#!/usr/bin/env python3
"""A short script to remove version details from src columns in s rows in GP file"""

# import os
# import sys
import fileinput
import re

def drop_version(line, sep="\t"):
    """Drop version info from src column."""
    cols_lst = line.split()
    cols_lst[1] = re.sub(r'\.\d+$', '', cols_lst[1])

    replaceline = sep.join(cols_lst)
    return replaceline


with fileinput.input() as infiles:
    for line in infiles:
            replace = drop_version(line)
            print(replace)
