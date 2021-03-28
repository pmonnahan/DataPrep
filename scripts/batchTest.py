#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: Perform homogeneity/batch testing to identify SNPs with inconsistent genotype frequencies across datasets.  Testing is performed within cases or controls and also within particular ancestries.  Ancestries must be inferred a priori via an ancestry inference software (e.g. RFMix) or available via self-reporting (e.g. race/ethnicity).
Requirements: 
 Takes the following arguments:
    -k: Required. Input key that has 6 labeled columns: FID IID  Dataset   Phenotype  Group.  FID and IID should match the names in the plink files, Dataset indicates the particular dataset and Group indicates the race/ethnicity.  Samples sizes will likely only be sufficient for testing of white, black, and/or latino
    -i: Required. Prefix to plink files (do not include the plink suffix')
    -o: output_directory
    -m: Minimum number of individuals in a dataset required for batch testing. default = 10
    -t: p-value threshold for batch effect test.  default = 0.005
    -C: chromosome to consider; default is to consider all chromosomes, but splitting calculation by chromosome can greatly increase speed
    -p: plink_path, default="plink"
    -r: string to be used to specify resources for plink. default='--memory 8000 --threads 1'
    -q: flag requesting that results are printed as they are generated.
    -c: threads
Note: Testing just the imputed case datasets for ADMIRAL for both white and black groups took 24 cores, 64gb, and ~48 run time for chromosome 1.
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
def calc_freqs(Data, dataset, pheno, group, plink_file, min_ind, outdir, plink, rsrc_str, chromosome):
    """For the current ancestry and phenotype group, calculate genotype frequencies for the specified dataset"""
    cdat = Data[Data.Dataset == dataset] #Data is already subsetted to include just individuals of a particular ancestry and phenotype
    fname = f"{outdir}/{dataset}.{group}.{pheno}.samples" #Label output file
    if cdat.shape[0] > min_ind: #Check if this dataset contains minimum number of individuals for inclusion in testing
        cdat[["FID", "IID"]].to_csv(fname, sep="\t", index=False) #Write sample names to file to be used in plink subsetting
        if chromosome != "all": #Calculate frequency for all chromosomes or just a single chromosome
            cmd = f"{plink} --bfile {plink_file} --chr {chromosome} --keep {fname} --hardy --maf 0.001 --out {outdir}/{dataset}.{group}.{pheno} {rsrc_str}"
        else:
            cmd = f"{plink} --bfile {plink_file} --keep {fname} --hardy --maf 0.001 --out {outdir}/{dataset}.{group}.{pheno} {rsrc_str}"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish
        if pheno == "case": cmd = f"grep -w AFF " #Make sure we grab the appropriate rows from PLINK output
        if pheno == "control": cmd = f"grep -w UNAFF "
        cmd += f"{outdir}/{dataset}.{group}.{pheno}.hwe | awk '{{print $2,$6}}' | sed 's:/:\t:g' > {outdir}/{dataset}.{group}.{pheno}.genofreqs"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish
    return()

def wrap_freq(Dat, plink_file, min_ind, outdir, plink, rsrc_str, threads, chromosome):
    """Loop over disease state and ethnicities to calculate genotype frequencies for each"""
    for pheno in ['case', 'control']:
        for group in set(Dat.Group): #Loop over groups present in input key
            if pheno == 'case': tDat = Dat[(Dat.Group == group) & (Dat.Phenotype == 2)] #Extract individuals corresponding to current group and phenotype
            else: tDat = Dat[(Dat.Group == group) & (Dat.Phenotype == 1)]
            if len(set(tDat.Dataset)) > 1: #Batch testing only makes sense when more than one dataset has individuals from current group and phenotype
                arg_list = []
                for dataset in set(tDat.Dataset): #Format argument list for each dataset
                    arg_list.append((tDat, dataset, pheno, group, plink_file, min_ind, outdir, plink, rsrc_str, chromosome))
                pool = multiprocessing.Pool(threads)
                pool.starmap(calc_freqs, arg_list) #Pass arguments to function that calculates genotype frequencies
                pool.close()
                pool.join()

                cmd = f"cat {outdir}/*{group}.{pheno}.genofreqs | sort -k 1 > {outdir}/{group}.{pheno}.genofreqs" #Concatenate genotype frequencies from same group and phenotype and Sort so that data from same SNP from different datasets are adjacent
                pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
                out1, err1 = pp1.communicate()  # Wait for it to finiss

    return(out1,err1)

def batchTest(snp_array, alternative, hybrid, midP_method, simulate_pval, replicates):
    """Perform fisher exact test on set of genotype frequencies for a specific SNP"""
    res = FisherExact.fisher_exact(snp_array, alternative, hybrid, midP_method, simulate_pval, replicates)
    return(res)

def wrapTest(genofreq_path, threads, names=["SNP", "AA", "AB", "BB"], sep=r"\s+", batchSize=100, reps = 2000, quickprint=False, threshold = 0.005):
    SNP = {} #dictionary to hold results of batchTest for each SNP
    arg_list = [] #Holds arguments for each SNP to be sent as batches to batchTest function
    ids = []
    batches = 0
    if quickprint: #Print results as we go?
        pfile = open(genofreq_path.replace("genofreqs","pvals"), 'w') #File that will hold all results
        efile = open(f"{genofreq_path}.exclude", 'w') #File that will hold SNPs whose p-value is below threshold.
    if os.path.exists(genofreq_path):
        Data = pd.read_csv(genofreq_path, names=names, sep=sep) #Read in concatenated genotype frequencies.
        for snp in set(Data.SNP): #Loop over SNPs in dataset
            dat = Data[Data.SNP == snp] #Extract data for this SNP
            snp_array = np.array(dat[["AA", "AB", "BB"]]) #Convert to numpy array
            # Catch and remove monomorphic SNPs
            if np.count_nonzero(snp_array[:,1]) > 0 and (np.count_nonzero(snp_array[:,2]) > 0 or np.count_nonzero(snp_array[:,0]) > 0):
                arg_list.append((snp_array,'two-sided', True, False, True, reps))
                ids.append(snp)
            if len(arg_list) % batchSize == 0: #Once a full batch of SNPs has been assembled into arg_list, submit for parallel processing and wait to finish
                pool = multiprocessing.Pool(threads)
                res = pool.starmap(batchTest, arg_list)
                pool.close()
                pool.join()
                if quickprint: #Print results of this batch
                    for var, pval in zip(ids, res):
                        pfile.write(f"{var}\t{pval}\n")
                        if pval < threshold: efile.write(f"{var}\n")
                else: #Or, hold all results in memory
                    for var, pval in zip(ids, res): SNP[var] = pval
                arg_list = [] #Reset arg_list to being storing arguments for next batch
                ids = []
                batches += 1
                print(f"{batches / ceil(len(set(Data.SNP))/batchSize)} done with {genofreq_path}")
        if arg_list: #Process any remaining SNPs if they did not constitute a full batch
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

Dat = pd.read_csv(args.k, sep = r"\s+") #Read in input key
min_ind = args.m
threshold = float(args.t)
plink_file = args.i
plink = args.p
rsrc_str = args.r
out = args.o
if args.q: quickprint = True
else: quickprint = False

if not os.path.exists(args.o): os.mkdir(args.o) #Create output directory if it does not exist

o1, e1 = wrap_freq(Dat, plink_file, min_ind, out, plink, rsrc_str, args.c, args.C) #Calculate genotype frequencies for each combination of disease status and ancestry

if not e1: #Ensure that genotype frequency calculation completed without error
    for pheno in ['case', 'ctrl']:
        for group in set(Dat.Group): #Loop over ancestries
        #for group in ["Black"]:
            genofreq_path = f"{out}/{group}.{pheno}.genofreqs"
            Results = wrapTest(genofreq_path, args.c, quickprint=quickprint, threshold=threshold) #Perform batch tests for this group/disease status
            if not quickprint: #Write results if they were not written as generated
                with open(f"{out}/{group}.{pheno}.pvals", 'w') as pfile:
                    with open(f"{out}/{group}.{pheno}.batchTest.exclude", 'w') as exclude:
                        for snp, pval in Results.items():
                            pfile.write(f"{snp}\t{pval}\n")
                            if pval < threshold: exclude.write(f"{snp}\n")
else:
    print("Calculation of genotype frequencies failed")




