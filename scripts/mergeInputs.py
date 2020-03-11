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
import itertools
import shutil
from functools import reduce
import multiprocessing


def flipstrand_and_updateID(input, names_file, snps2flip, output, plink):
    cmd = f"{plink} --bfile {input} --update-name {names_file} --flip {snps2flip} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()
    return(output)


def var_flt(bfile, keptVariants, output, plink):
    cmd = f"{plink} --bfile {bfile} --extract {keptVariants} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(output)


def parallelQC(arg_list, threads, function="flipstrand_and_updateID"):
    pool = multiprocessing.Pool(threads)
    if function == "flipstrand_and_updateID":
        outfiles = pool.starmap(flipstrand_and_updateID, arg_list)
    if function == "var_flt":
        outfiles = pool.starmap(var_flt, arg_list)
    pool.close()
    pool.join()
    return(outfiles)


def merge_files(file, file_number, output, samples='all', plink='plink'):
    if file_number == 1:
        if samples == 'all': cmd = f"{plink} --bfile {file} --make-bed --allow-no-sex --keep-allele-order --out {output}"
        else: cmd = f"{plink} --bfile {file} --keep {samples} --make-bed --allow-no-sex --out {output}"
    else:
        if samples == 'all': cmd = f"{plink} --merge-list {file} --make-bed --allow-no-sex --out {output}"
        else: cmd = f"{plink} --merge-list {file} --keep {samples} --make-bed --allow-no-sex --out {output}"
    print(cmd)
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)


def read_bims(file_key):
    DF_dict = {}
    with open(file_key, 'r') as files:
        for file in files:
            df = pd.read_csv(f"{file.strip()}.bim", sep=r"\s+", header=None, dtype=str)
            DF_dict[os.path.basename(file.strip())] = df
    return(DF_dict)


def flip_scan(input, output, plink):
    cmd = f"{plink} --bfile {input} --flip-scan --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()

    return (output)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma separated list(no spaces) with the dataset PREFIXes.  Datasets should be in plink format.  ')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-t', type=int, metavar='threads', default=1, help='')
    parser.add_argument('-s', type=str, metavar='sample_file', default='all',
                        help='File containing sample IDs to retain from merged plink files. One sample per line')
    args = parser.parse_args()

    files = args.i.split(",")

    indir_dict = {}
    with open(f"{args.d}/files2merge.txt", 'w') as merge_file:
        for file in files:
            merge_file.write(f"{args.d}/{file}\n")
            indir_dict[os.path.basename(file.strip())] = os.path.dirname(file)

    DF_dict = read_bims(f"{args.d}/files2merge.txt")

    # Identify merge issues before they occur.  SNPs with matching positions but different names.  Matching SNPs where strand is flipped. e.g. A-G in one is T-C in the other
    Names_dict = {}
    Flips_dict = {}
    for i, comb in enumerate(itertools.combinations(DF_dict.keys(), r=2)):
        if i==0:  # Hold first data frame as the reference and flip needed sites in subsequent data frames with respect to constant reference
            reference = comb[0]
        df1 = DF_dict[comb[0]]
        df2 = DF_dict[comb[1]]
        df = pd.merge(df1, df2, on=(0,3))
        pdb.set_trace()
        # TODO: need to catch when names match, alleles match, but positions do not.
        bad_names = df[df['1_x'] != df['1_y']]
        bad_names_flips = bad_names[(bad_names['4_x'] != bad_names['4_y']) & (bad_names['4_x'] != bad_names['5_y']) & (bad_names['5_x'] != bad_names['4_y']) & (bad_names['5_x'] != bad_names['5_y'])]

        df_b = pd.merge(df1, df2, on=1)
        good_names_flips = df_b[(df_b['4_x'] != df_b['4_y']) & (df_b['4_x'] != df_b['5_y']) & (df_b['5_x'] != df_b['4_y']) & (df_b['5_x'] != df_b['5_y'])]

        if comb[0] == reference:
            Names_dict[comb[1]] = bad_names[['1_y', '1_x']]
            flips = list(good_names_flips[1]) + list(bad_names_flips['1_x'])  # SNP ids to be flipped should be specified in terms of reference ids bc
            Flips_dict[comb[1]] = flips

    # harmonize SNPids across data sets
    with open(f"{args.d}/files2premerge.txt", 'w') as mergefile:
        mergefile.write(f"{indir_dict[reference]}/{reference}\n")
        arg_list4 = []
        for name, data in Names_dict.items():
            data[['1_y', '1_x']].to_csv(f"{args.d}/{name}_id_updates", sep="\t", index=False, header=False)
            with open(f"{args.d}/{name}_flips", 'w') as flip_file:
                flip_file.write('\n'.join(Flips_dict[name]))
            arg_list4.append((f"{indir_dict[name]}/{name}", f"{args.d}/{name}_id_updates", f"{args.d}/{name}_flips", f"{args.d}/{name}.ids", args.p))
            mergefile.write(f"{args.d}/{name}.ids\n")

    outs = parallelQC(arg_list4, args.t, function="flipstrand_and_updateID")
    # pdb.set_trace()
    DF_dict = read_bims(f"{args.d}/files2premerge.txt")

    #Losing LOTS of SNPs here.
    # Poor overlap between aric and stjude/cog9906
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on=1), DF_dict.values())
    df[1].to_csv(f"{args.d}/mergeSNPs.txt", header=False, index=False)

    arg_list5 = []
    with open(f"{args.d}/files2merge.txt", 'w') as mergefile:
        with open(f"{args.d}/files2premerge.txt", 'r') as premerge:
            for input in premerge:
                input = input.strip()
                arg_list5.append((input, f"{args.d}/mergeSNPs.txt", f"{input}.flt", args.p))
                mergefile.write(f"{input}.flt\n")

    outs = parallelQC(arg_list5, args.t, function="var_flt")

    # First merge attempt is likely to have strand flip errors
    out1, err1 = merge_files(f"{args.d}/files2merge.txt", len(files), f"{args.d}/{args.o}", samples=args.s, plink=args.p)

    # Restore Ref/Alt alleles using bim file of first dataset?



    # TODO: perform flip-scan via plink to look for incorrect strand assignment in a subset of sample after finding best-flipped combon








