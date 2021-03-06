#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: 
Requirements: 
 Takes the following arguments:
    -i: 
    -o:
Date:
"""

# Import Modules
import allel
import subprocess
import argparse
import os
import pdb

# Define Functions


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()

    with open(args.i, 'r') as input_file:
        for line in input_file:
            if not line.startswith("#"):
                line.replace('"','').strip(",")
                if not line[0] == "Probe Set ID":
                    if line[4] == "+"




