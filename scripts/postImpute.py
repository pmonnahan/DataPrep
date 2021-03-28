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
import subprocess
import argparse
import os
import pdb
import pandas as pd

# Define Functions
def badVars(infile, threshold):
    dat = pd.read_csv(infile, compression='gzip', sep = r"\s+")
    ndat = dat.query(f"Rsq < {threshold}")
    return(ndat.SNP)


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma-separated string (no space) of the directories containing the info.gz files that are output alongsize the dosage.vcf.gz files')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-r', type=str, metavar='treshold_Rsquared', default=0.7, help='imputed variants with Rsquared values lower than this threshold will be removed')
    parser.add_argument('-g', type=str, metavar='treshold_GenotypeProbs', default=0.85,
                        help='Set genotypes below this to change to missing')
    parser.add_argument('-P', type=str, metavar='file prefix', default="chr",
                        help='prefix preceeding the chromosome number')
    parser.add_argument('-c', type=str, metavar='chromosome_number', default="all",
                        help='set to \'all\' to loop over autosomes')

    args = parser.parse_args()

    in_dirs = args.i.split(",")

    if args.c == 'all':
        chroms = [x for x in range(1, 23)]
    else: chroms = [args.c]

    for chrom in chroms:
        for i, in_dir in enumerate(in_dirs):
            file = f"{in_dir}/{args.P}{chrom}.info.gz"
            if os.path.exists(file):
                bad_vars = badVars(file, args.r)
                if i == 0:
                    DF = bad_vars
                else:
                    DF = pd.concat([DF, bad_vars])
            else:
                print(f"Did not find {file}")
        DF.drop_duplicates(inplace=True)
        DF.to_csv(f"{args.d}/chr{chrom}_fltSNPs.txt", header=False, index=False)
        with open(f"{args.d}/2bMerged.txt", 'w') as file_list:
            for i, in_dir in enumerate(in_dirs):
                vcf = f"{in_dir}/{args.P}{chrom}.dose.vcf.gz"
                cmd = f"{args.p} --vcf {vcf} --vcf-min-gp {args.g} --exclude {args.d}/chr{chrom}_fltSNPs.txt --keep-allele-order --make-bed --out {in_dir}/chr{chrom}"
                pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
                out1, err1 = pp1.communicate()  # Wait for it to finish
                file_list.write(f"{in_dir}/chr{chrom}\n")

        cmd = f"{args.p} --bmerge {in_dir}/chr{chrom} --keep-allele-order --make-bed --out {args.d}/{args.o}"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish







