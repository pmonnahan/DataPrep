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
from distutils import spawn
import multiprocessing

# Define Functions
def var_miss(bfile, threshold, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --geno {threshold} --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def geno_miss(bfile, threshold, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --mind {threshold} --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def geno_flt(bfile, keptSamples, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --keep {keptSamples} --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def var_flt(bfile, keptVariants, output, plink, rsrc_str):
    cmd = f"{plink} --bfile {bfile} --extract {keptVariants} --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def calc_stats(bfile, plink, rsrc_str):
    '''calculate frequency, missingness, HWE prob, and missingness by sex'''

    cmd = f"{plink} --bfile {bfile} --freq --missing --hardy --mpheno 3 --pheno {bfile}.fam --test-missing --out {bfile} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish

    return(out1, err1)

def wrapQC(input, output, tvm1, tgm, tvm2, maf, hwe, mbs, plink, rsrc_str, snps_only=False):

    out1, err1 = var_miss(input, tvm1, f"{output}.var_miss{tvm1}", plink, rsrc_str)
    out2, err2 = geno_miss(f"{output}.var_miss{tvm1}", tgm, f"{output}.geno_miss{tgm}", plink, rsrc_str)
    out3, err3 = geno_flt(input, f"{output}.geno_miss{tgm}.fam",
                         f"{output}.geno_flt{tgm}", plink, rsrc_str)

    out2, err2 = calc_stats(f"{output}.geno_flt{tgm}", plink, rsrc_str)
    DF_list = []
    flt_strings = {'frq': f"frq_MAF < {maf}", 'hwe': f"hwe_P < {hwe}", 'lmiss': f"lmiss_F_MISS > {tvm2}", 'missing': f"missing_P < {mbs}"}
    for stat in ['frq', 'hwe', 'lmiss', 'missing']:
        try:
            tdf = pd.read_csv(f"{output}.geno_flt{tgm}.{stat}", sep = r"\s+")
            tdf = tdf.add_prefix(f"{stat}_")
            tdf['SNP'] = tdf[f"{stat}_SNP"]
            if stat=="hwe":
                tdf = tdf[tdf['hwe_TEST']=="UNAFF"]
            DF_list.append(tdf)
            tdf = tdf.query(flt_strings[stat])
            tdf[f"{stat}_SNP"].to_csv(f"{output}.fltdSNPs.{stat}.txt", header=False, index=False)
            tdf['SNP'].to_csv(f"{output}.fltdSNPs.txt", header=False, index=False, mode='a')
        except FileNotFoundError:
            print(f"Did not find file: {output}.geno_flt{tgm}.{stat}")
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on='SNP'), DF_list)
    if snps_only:
        awk_string = "awk '{split(\"A T C G\",tmp); for (i in tmp) arr[tmp[i]]; if(!($5 in arr) || !($6 in arr)) print $2}' >> "
    else:
        awk_string = "awk '{if($5 == \"N\" || $6 == \"N\") print $2}' >> "
    cmd = f"{plink} --bfile {output}.geno_flt{tgm} --list-duplicate-vars suppress-first --out {output} {rsrc_str} ; tail -n +2 {output}.dupvar | cut -f 4 >> {output}.fltdSNPs.txt;"
    cmd += f"cat {output}.geno_flt{tgm}.bim | " + awk_string + output + ".fltdSNPs.txt; "
    cmd += f"cut -f 2 {output}.geno_flt{tgm}.bim | sort | uniq -d >> {output}.fltdSNPs.txt; {plink} --bfile {output}.geno_flt{tgm} --exclude {output}.fltdSNPs.txt --make-bed --out {output} {rsrc_str}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish

    return(output)


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma-separated string (no space) of the PLINK files to be QCed and combined')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-tvm1', type=str, metavar='treshold_variant_missingness1', default=0.2, help='Used for generated filtered dataset for calculating genotype missingness')
    parser.add_argument('-tgm', type=str, metavar='treshold_genotype_missingness', default=0.1, help='')
    parser.add_argument('-tvm2', type=str, metavar='treshold_variant_missingness2', default=0.05, help='Missingness allowed for variants post genotype filtration')
    parser.add_argument('-maf', type=str, metavar='minor_allele_frequency', default=0.01,
                        help='')
    parser.add_argument('-mbs', type=str, metavar='missingness_by_sex', default=0.0000001,
                        help='')
    parser.add_argument('-hwe', type=str, metavar='p_value_cutoff_HWEtest', default=-1,
                        help='')
    parser.add_argument('-r', type=str, metavar='resource_string', default='1',
                        help='string to be used to specify resources for plink. e.g. --memory 16000 --threads 1'),
    parser.add_argument('-s', type=str, metavar='sample_file', default='all',
                        help='File containing sample IDs to retain from merged plink files. One sample per line')
    parser.add_argument('--snps_only', action='store_true', help="output only snps; i.e. filter indels")
    args = parser.parse_args()


    if not os.path.exists(args.d):
        os.mkdir(args.d)

    if args.p == "plink":
        args.p = spawn.find_executable(args.p)
    if args.snps_only:
        wrapQC(args.i, f"{args.d}/{args.o}", args.tvm1, args.tgm, args.tvm2, args.maf, args.hwe, args.mbs, args.p, args.r, snps_only=True)
    else:
        wrapQC(args.i, f"{args.d}/{args.o}", args.tvm1, args.tgm, args.tvm2, args.maf, args.hwe, args.mbs, args.p, args.r)