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


def flipstrand_and_updateID(input, names_file, snps2flip, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {input} --update-name {names_file} --flip {snps2flip} --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()
    return(output)


def var_flt(bfile, keptVariants, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --extract {keptVariants} --keep-allele-order --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(output)


def exclude_variants(bfile, excludedVariants, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --exclude {excludedVariants} --keep-allele-order --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(output)


def parallelQC(arg_list, threads, function="flipstrand_and_updateID"):
    pool = multiprocessing.Pool(threads)
    if function == "flipstrand_and_updateID":
        outfiles = pool.starmap(flipstrand_and_updateID, arg_list)
    if function == "var_flt":
        outfiles = pool.starmap(var_flt, arg_list)
    if function == "exclude_variants":
        outfiles = pool.starmap(exclude_variants, arg_list)
    pool.close()
    pool.join()
    return(outfiles)


def merge_files(file, file_number, output, genotype_missing, samples='all', plink='plink', rsrc_str = '--memory 16000 --threads 1'):
    if file_number == 1:
        if samples == 'all': cmd = f"{plink} --bfile {file} --mind {genotype_missing} --make-bed --allow-no-sex --keep-allele-order --out {output} {rsrc_str}"
        else: cmd = f"{plink} --bfile {file} --keep {samples} --mind {genotype_missing} --keep-allele-order --make-bed --allow-no-sex --out {output} {rsrc_str}"
    else:
        if samples == 'all': cmd = f"{plink} --merge-list {file} --mind {genotype_missing} --keep-allele-order --make-bed --allow-no-sex --out {output} {rsrc_str}"
        else: cmd = f"{plink} --merge-list {file} --keep {samples} --mind {genotype_missing} --keep-allele-order --make-bed --allow-no-sex --out {output} {rsrc_str}"
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


def flip_scan(input, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {input} --flip-scan --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()

    return (output)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge together a set of plink files.  Assumes that reference allele has been harmonized across all datasets.')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma separated list(no spaces) with the dataset PREFIXes.  Datasets should be in plink format.  ')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-m', type=str, metavar='genotype_missingness', default="0.1", help='')
    parser.add_argument('-s', type=str, metavar='sample_file', default='all',
                        help='File containing sample IDs to retain from merged plink files. One sample per line')
    parser.add_argument('-r', type=str, metavar='resource_string', default='--memory 16000 --threads 1',
                        help='string to be used to specify resources for plink. e.g. --memory 16000 --threads 1'),
    parser.add_argument('--drop_pos_dups', action='store_true', help="remove any variants with positional duplicates. Designed to catch multiallelic variants that have been specified as multiline biallelic variants")
    args = parser.parse_args()

    files = args.i.split(",")

    indir_dict = {}
    with open(f"{args.d}/{args.o}_files2merge.txt", 'w') as merge_file: # Do dup dropping as optional step in here and then can sub modded file into merge file spot
        for file in files:
            if args.drop_pos_dups:
                df = pd.read_csv(f"{file.strip()}.bim", sep=r"\s+", header=None, dtype=str)
                df.drop_duplicates(3, keep=False, inplace=True)
                df[1].to_csv(f"{args.d}/{file.strip()}_unique", sep="\t", index=False, header=False)
                var_flt(f"{args.d}/{file.strip()}", f"{args.d}/{file.strip()}_unique", f"{args.d}/{file.strip()}.unq", args.p, args.r)
                merge_file.write(f"{args.d}/{file.strip()}.unq\n")
                indir_dict[os.path.basename(file.strip()) + ".unq"] = os.path.dirname(file)
            else:
                merge_file.write(f"{args.d}/{file.strip()}\n")
                indir_dict[os.path.basename(file.strip())] = os.path.dirname(file)

    DF_dict = read_bims(f"{args.d}/{args.o}_files2merge.txt")

    # Identify merge issues before they occur.  SNPs with matching positions but different names.  Matching SNPs where strand is flipped. e.g. A-G in one is T-C in the other
    Names_dict = {}
    Flips_dict = {}
    Pos_dict = {}  # If dropping positional duplicates, then updating variant IDs and strand is not necessary.  HOWEVER, the attempts above, to identify and remove positional duplicates resulting from multiline multiallelic variants will miss instances where the alternative allele differes across populations.  That is, the variant is biallelic within each population, but with a different alternative allele.  In these cases, we want to identify all naming conflicts and remove them.
    for i, comb in enumerate(itertools.combinations(DF_dict.keys(), r=2)):
        if i==0:  # Hold first data frame as the reference and flip needed sites in subsequent data frames with respect to constant reference
            reference = comb[0]
        # pdb.set_trace()
        df1 = DF_dict[comb[0]]
        df2 = DF_dict[comb[1]]
        df = pd.merge(df1, df2, on=(0, 3))  # Merging on chromosome and position
        # pdb.set_trace()
        # TODO: need to catch when names match, alleles match, but positions do not.
        bad_names = df[df['1_x'] != df['1_y']]  # Check if variant IDs match
        # Check if, for variants with matching position and non-matching names, are ALL alleles the different (indicating a flip strand).
        bad_names_flips = bad_names[(bad_names['4_x'] != bad_names['4_y']) & (bad_names['4_x'] != bad_names['5_y']) & (bad_names['5_x'] != bad_names['4_y']) & (bad_names['5_x'] != bad_names['5_y'])]

        df_b = pd.merge(df1, df2, on=1)
        good_names_flips = df_b[(df_b['4_x'] != df_b['4_y']) & (df_b['4_x'] != df_b['5_y']) & (df_b['5_x'] != df_b['4_y']) & (df_b['5_x'] != df_b['5_y'])]
        good_names_all_conflicts = df_b[(df_b['4_x'] != df_b['4_y']) | (df_b['5_x'] != df_b['5_y'])]
        # Record all name conflicts in case we want to delete all positional duplicates later
        try:
            Pos_dict[comb[1]] += list(bad_names['1_y'])
            Pos_dict[comb[0]] += list(bad_names['1_x'])
            Pos_dict[comb[1]] += list(good_names_all_conflicts[1])
            Pos_dict[comb[0]] += list(good_names_all_conflicts[1])
        except (KeyError, AttributeError):
            Pos_dict[comb[1]] = list(bad_names['1_y'])
            Pos_dict[comb[0]] = list(bad_names['1_x'])
            Pos_dict[comb[1]] += list(good_names_all_conflicts[1])
            Pos_dict[comb[0]] += list(good_names_all_conflicts[1])

        if comb[0] == reference:
            Names_dict[comb[1]] = bad_names[['1_y', '1_x']]
            flips = list(good_names_flips[1]) + list(bad_names_flips['1_x'])  # SNP ids to be flipped should be specified in terms of reference ids bc
            Flips_dict[comb[1]] = flips

    # harmonize SNPids across data sets OR just exclude name conflicts (which will be multiline multiallelic variants)
    with open(f"{args.d}/{args.o}_files2premerge.txt", 'w') as mergefile:
        arg_list4 = []
        # pdb.set_trace()
        if not args.drop_pos_dups:
            mergefile.write(f"{indir_dict[reference]}/{reference}\n")
            for name, data in Names_dict.items():

                data[['1_y', '1_x']].to_csv(f"{indir_dict[name]}/{name}_id_updates", sep="\t", index=False, header=False)
                with open(f"{indir_dict[name]}/{name}_flips", 'w') as flip_file:
                    flip_file.write('\n'.join(Flips_dict[name]))
                arg_list4.append((f"{indir_dict[name]}/{name}", f"{indir_dict[name]}/{name}_id_updates", f"{indir_dict[name]}/{name}_flips", f"{indir_dict[name]}/{name}.ids", args.p, args.r))
                mergefile.write(f"{indir_dict[name]}/{name}.ids\n")
        else:
            for name, data in Pos_dict.items():
                if Pos_dict[name]:
                    with open(f"{indir_dict[name]}/{name}_dup", 'w') as dup_file:
                        dup_file.write('\n'.join(Pos_dict[name]))
                    arg_list4.append((f"{indir_dict[name]}/{name}", f"{indir_dict[name]}/{name}_dup", f"{indir_dict[name]}/{name}.unq2", args.p, args.r))
                    mergefile.write(f"{indir_dict[name]}/{name}.unq2\n")
                else: mergefile.write(f"{indir_dict[name]}/{name}\n")

    if args.drop_pos_dups and arg_list4:
        outs = parallelQC(arg_list4, 1, function="exclude_variants")
    else:
        outs = parallelQC(arg_list4, 1, function="flipstrand_and_updateID")
    # pdb.set_trace()
    DF_dict = read_bims(f"{args.d}/{args.o}_files2premerge.txt")

    #Losing LOTS of SNPs here.
    # Poor overlap between aric and stjude/cog9906
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on=1), DF_dict.values())
    df[1].to_csv(f"{args.d}/{args.o}_mergeSNPs.txt", header=False, index=False)

    arg_list5 = []
    with open(f"{args.d}/{args.o}_files2merge.txt", 'w') as mergefile:
        with open(f"{args.d}/{args.o}_files2premerge.txt", 'r') as premerge:
            for input in premerge:
                input = input.strip()
                arg_list5.append((input, f"{args.d}/{args.o}_mergeSNPs.txt", f"{input}.flt", args.p, args.r))
                mergefile.write(f"{input}.flt\n")

    outs = parallelQC(arg_list5, 1, function="var_flt")

    # First merge attempt is likely to have strand flip errors
    out1, err1 = merge_files(f"{args.d}/{args.o}_files2merge.txt", len(files), f"{args.d}/{args.o}", args.m, samples=args.s, plink=args.p, rsrc_str=args.r)

    # Restore Ref/Alt alleles using bim file of first dataset?



    # TODO: perform flip-scan via plink to look for incorrect strand assignment in a subset of sample after finding best-flipped combon








