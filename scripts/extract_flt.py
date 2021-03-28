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
# import allel
import subprocess
# import os
# import pdb
import pandas as pd

# Define Functions
mbc_threshold = 0.000001
hwe_threshold = 0.000001
mis_threshold = 0.05
input = "ADMIRAL6-QC"
euro_samps = "/home/spectorl/pmonnaha/misc/EUR_BALL_CaseControl.txt"
cmd = f"module load plink; plink --bfile {input} --keep {euro_samps} --keep-allele-order --make-bed --out {input}.EUR"
pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
out1, err1 = pp1.communicate()  # Wait for it to finish
cmd = f"module load plink; plink --bfile {input}.EUR --test-missing --allow-no-sex --out {input}.EUR"
pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
out1, err1 = pp1.communicate()  # Wait for it to finish
cmd = f"module load plink; plink --bfile {input}.EUR --missing --hardy --allow-no-sex --out {input}.EUR"
pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
out1, err1 = pp1.communicate()  # Wait for it to finish
mbc = pd.read_table(f"{input}.EUR.missing", sep=r"\s+")
mis = pd.read_table(f"{input}.EUR.lmiss", sep=r"\s+")
hwe = pd.read_table(f"{input}.EUR.hwe", sep=r"\s+")
hwe = hwe[hwe['TEST']=="UNAFF"]
dat = pd.concat([mbc[mbc['P'] < mbc_threshold]['SNP'],
                 mis[mis['F_MISS'] > mis_threshold]['SNP'], hwe[hwe['P'] < hwe_threshold]['SNP']]).drop_duplicates()
dat.to_csv(f'{input}.EUR.filter.txt', header=False, index=False)

cmd = f"module load plink; plink --bfile {input}.EUR --exclude {input}.EUR.filter.txt --keep-allele-order --make-bed --out {input}.EUR.QC"
pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
out1, err1 = pp1.communicate()  # Wait for it to finish

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()

