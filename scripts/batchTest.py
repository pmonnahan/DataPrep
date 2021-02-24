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
import numpy as np
import multiprocessing
import FisherExact
from math import ceil

# Define Functions
def calc_freqs(Cases, dataset, group, plink_file, min_ind, outdir, plink, rsrc_str, chromosome):
    cdat = Cases[Cases.Dataset == dataset]
    fname = f"{outdir}/{dataset}.{group}.samples"
    if cdat.shape[0] > min_ind:
        cdat[["FID", "IID"]].to_csv(fname, sep="\t", index=False)
        if chromosome != "all":
            cmd = f"{plink} --bfile {plink_file} --chr {chromosome} --keep {fname} --hardy --maf 0.001 --out {outdir}/{dataset}.{group}.case {rsrc_str}"
        else:
            cmd = f"{plink} --bfile {plink_file} --keep {fname} --hardy --maf 0.001 --out {outdir}/{dataset}.{group}.case {rsrc_str}"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish
        cmd = f"grep ALL {outdir}/{dataset}.{group}.case.hwe | awk '{{print $2,$6}}' | sed 's:/:\t:g' > {outdir}/{dataset}.{group}.case.genofreqs"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish
    return()

def wrap_freq(Dat, plink_file, min_ind, outdir, plink, rsrc_str, threads, chromosome):
    for pheno in ['case','control']:
        for group in set(Dat.Group):
            if pheno == 'case': tDat = Dat[(Dat.Group == group) & (Dat.Phenotype == 2)]
            else: tDat = Dat[(Dat.Group == group) & (Dat.Phenotype == 1)]
            if len(set(tDat.Dataset)) > 1:
                arg_list = []
                for dataset in set(tDat.Dataset):
                    arg_list.append((tDat, dataset, group, plink_file, min_ind, outdir, plink, rsrc_str, chromosome))
                pool = multiprocessing.Pool(threads)
                pool.starmap(calc_freqs, arg_list)
                pool.close()
                pool.join()

                cmd = f"cat {outdir}/*{group}.{pheno}.genofreqs | sort -k 1 > {outdir}/{group}.{pheno}.genofreqs"
                pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
                out1, err1 = pp1.communicate()  # Wait for it to finiss

    return(out1,err1)

def batchTest(snp_array, alternative, hybrid, midP_method, simulate_pval, replicates):
    res = FisherExact.fisher_exact(snp_array, alternative, hybrid, midP_method, simulate_pval, replicates)
    return(res)

def wrapTest(genofreq_path, threads, names=["SNP", "AA", "AB", "BB"], sep=r"\s+", batchSize=100, reps = 2000, quickprint=False, threshold = 0.005):
    SNP = {}
    arg_list = []
    ids = []
    batches = 0
    if quickprint:
        pfile = open(genofreq_path.replace("genofreqs","pvals"), 'w')
        efile = open(f"{genofreq_path}.exclude", 'w')
    if os.path.exists(genofreq_path):
        Data = pd.read_csv(genofreq_path, names=names, sep=sep)
        for snp in set(Data.SNP):
            dat = Data[Data.SNP == snp]
            snp_array = np.array(dat[["AA", "AB", "BB"]])
            # Catch and remove monomorphic SNPs
            if np.count_nonzero(snp_array[:,1]) > 0 and (np.count_nonzero(snp_array[:,2]) > 0 or np.count_nonzero(snp_array[:,0]) > 0):
                arg_list.append((snp_array,'two-sided', True, False, True, reps))
                ids.append(snp)
            if len(arg_list) % batchSize == 0:
                pool = multiprocessing.Pool(threads)
                res = pool.starmap(batchTest, arg_list)
                pool.close()
                pool.join()
                if quickprint:
                    for var, pval in zip(ids, res):
                        pfile.write(f"{var}\t{pval}\n")
                        if pval < threshold: efile.write(f"{var}\n")
                else:
                    for var, pval in zip(ids, res): SNP[var] = pval
                arg_list = []
                ids = []
                batches += 1
                print(f"{batches / ceil(len(set(Data.SNP))/batchSize)} done with {genofreq_path}")
        if arg_list:
            pool = multiprocessing.Pool(threads)
            res = pool.starmap(batchTest, arg_list)
            pool.close()
            pool.join()
            if quickprint:
                for var, pval in zip(ids, res):
                    pfile.write(f"{var}\t{pval}\n")
                    if pval < threshold: efile.write(f"{var}\n")
                pfile.close()
                efile.close()
            else:
                for var, pval in zip(ids, res): SNP[var] = pval
    return(SNP)


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-k', type=str, metavar='key', required=True, help='input key that has 6 labeled columns: FID IID  Dataset   Phenotype  Group.  FID and IID should match the names in the plink files, Dataset indicates the particular dataset and Group indicates the race/ethnicity.  Samples sizes will likely only be sufficient for testing of white, black, and/or latino')
    parser.add_argument('-i', type=str, metavar='plink_file', required=True, help='prefix to plink files (do not include the plink suffix')
    parser.add_argument('-m', type=int, metavar='min_individuals', default=10, help='Minimum number of individuals in a dataset required for batch testing')
    parser.add_argument('-t', type=str, metavar='pval_threshold', default="0.005", help='p-value threshold for batch effect test')
    parser.add_argument('-C', type=str, metavar='chromosome', default="all",
                        help='chromosome to consider; default is to consider all chromosomes, but splitting calculation by chromosome can greatly increase speed')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-r', type=str, metavar='resource_string', default='--memory 8000 --threads 1',
                        help='string to be used to specify resources for plink. e.g. --memory 8000 --threads 1')
    parser.add_argument('-o', type=str, metavar='output_directory', default="", help='')
    parser.add_argument('-q', action="store_true", help='print results as they are generated')
    parser.add_argument('-c', type=int, metavar='threads', default=1, help='')
    args = parser.parse_args()

Dat = pd.read_csv(args.k, sep = r"\s+")
min_ind = args.m
threshold = float(args.t)
plink_file = args.i
plink = args.p
rsrc_str = args.r
out = args.o
if args.q: quickprint = True
else: quickprint = False

if not os.path.exists(args.o): os.mkdir(args.o)

o1, e1 = wrap_freq(Dat, plink_file, min_ind, out, plink, rsrc_str, args.c, args.C)

#TODO write separate function (or a quick print option) that prints to file as tests complete.  No reason to store the results the whole time.

if not e1:
    for pheno in ['case', 'ctrl']:
        for group in set(Dat.Group):
        #for group in ["Black"]:
            genofreq_path = f"{out}/{group}.{pheno}.genofreqs"
            Results = wrapTest(genofreq_path, args.c, quickprint=quickprint, threshold=threshold)
            if not quickprint:
                with open(f"{out}/{group}.{pheno}.pvals", 'w') as pfile:
                    with open(f"{out}/{group}.{pheno}.batchTest.exclude", 'w') as exclude:
                        for snp, pval in Results.items():
                            pfile.write(f"{snp}\t{pval}\n")
                            if pval < threshold: exclude.write(f"{snp}\n")
else:
    print("Calculation of genotype frequencies failed")




