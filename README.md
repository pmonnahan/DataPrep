# Pre-imputation QC pipeline

The purpose of this pipeline is to perform the following for a set of input PLINK datasets:

-   basic QC (genotype/variant missingness, HWE, and minor allele frequency)
-   harmonize allele specifications with the GRCh37 reference genome
-   produce a set of VCF files (separated by chromosome) for imputation
-   merge filtered datasets into a single dataset consisting only of overlapping sites.

A companion pipeline, which performs post-imputation QC, will download alongside the pre-imputation pipeline.  To use the post-imputation pipeline, see the README in the postImpute directory.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

 - [Requirements](#requirements)
   - [Snakemake](#snakemake)
   - [Singularity](#singularity)
 - [Running the workflow](#running-the-workflow)
   - [Other Notes](#other-notes)
    - [Debugging and error reports](#debugging-and-error-reports)
 - [Pipeline Overview](#pipeline-overview)
   - [Input Data](#input-data)
   - [Output](#output)
   - [Data Harmonization](#dataset-harmonization)
   - [Reference allele fixing](#reference-allele-fixing)
   - [Basic QC](#basic-qc)
   - [Merging Inputs (Optional)](#merging-inputs-optional)
   - [Imputation Preparation](#imputaton-preparation)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

![Pipeline DAG](https://github.com/pmonnahan/DataPrep/blob/master/Pipeline_DAG.png)

## Requirements

### Snakemake
The pipeline is coordinated and run on an HPC (or locally) using _Snakemake_.  On UMN's HPC, snakemake can be installed by:

    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge -c bioconda snakemake python=3.6

The 'module load' command will likely need to be run each time prior to use of Snakemake.

Alternatively, you can try installing snakemake via _pip_:

    pip3 install --user snakemake pyaml

### Singularity

The installation of the individual programs  used throughout this pipeline can be completely avoid by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image [here](https://github.com/pmonnahan/AncInf/blob/master/singularity/Singularity_defs.def).  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

However, in order to utilize the singularity image, _Singularity_ must be installed on the HPC.  Currently, the pipeline assumes that _Singularity_ will be available as a module and can be loaded into the environment via the command specified in the config.yml file, where it says 'singularity_module'.  The default setting will work for MSI at UMN.

Singularity settings in config.yml

    singularity:
      use_singularity: 'true'
      image: '/home/pmonnaha/pmonnaha/singularity/AncestryInference.sif

## Running the workflow

Clone the parent repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/DataPrep.git preImputeQC
    
The critical files responsible for executing the pipeline are contained in the *./workflow* subdirectory contained within the cloned repo.  They are: 

* Snakefile
* config.yml
* cluster.yaml  

The **Snakefile** is the primary workhouse of _snakemake_, which specifies the dependencies of various parts of the pipeline and coordinates execution.  No modifications to the **Snakefile** are necessary.  

In order for the **Snakefile** to locate all of the necessary input and correctly submit jobs to the cluster, **both** the **config.yaml** and **cluster.yaml** need to be modified. Open these files and change the required entries that are indicated with 'MODIFY'.  Other fields do not require modification, although this may be desired given the particulars of the run you wish to implement.  Details on each entry in the config file (e.g. what the program expects in each entry as well as the purpose of the entry) are provided in the _Pipeline Overview_ at the bottom.

The entire pipeline can be executed on a local machine (not recommended) or on an HPC, and the **cluster.yaml** file is required only for the latter.  For a local run, change the 'local_run' entry to 'true' under the 'run_settings' section of the config file, and launch snakemake from within the parent directory by the simple command:

    snakemake

However, multiple steps in the pipeline have high resource demands, and so are unlikely to be able to be run locally.  This option exists primarily for testing and troubleshooting, so the remainder of the  documentation assumes that the pipeline will be executed on an HPC.  In order to coordinate the use of the HPC, the following modifications to the snakemake command are required:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32

where -j specifies the number of jobs that can be submitted at once.  Note that the 'qsub' command is specific to the commonly-used PBS scheduler.  To run on a different HPC scheduler, the command would need to be modified accordingly.  For example, to coordinate submission to a slurm scheduler, the following command would be used:

    snakemake --cluster "sbatch --no-requeue --time={cluster.time} --mem-per-cpu={cluster.mem-per-cpu} --ntasks={cluster.ntasks} --nodes={cluster.nodes} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -o {cluster.o} -e {cluster.e} -A {cluster.A}"" --cluster-config workflow/cluster_yale.yaml -j 32

Note also that a different **cluster.yaml** file is required for the different scheduler.  If you open and inspect the **cluster.yaml** file vs the **cluster_yale.yaml** file, you will see syntax that is specific to PBS and slurm schedulers, respectively.  

One additional change in the **config.yml** is needed in order to correctly submit jobs to the HPC.  The relevant entries are under the 'run_settings' section of the config file:

    run_settings:
      local_run: 'false'
      cluster_config: 'workflow/cluster.yaml'
      scheduler: 'pbs'
      
Here, it is necessary that the 'cluster_config' entry is set to the path of the cluster.yaml file that will be used in the snakemake command.  Also, the scheduler must correspond to the syntax used in the snakemake command and cluster.yaml file.  I should point out that these additional changes are needed for responsibly using PLINK within a snakemake framework, and are not directly needed for snakemake.  PLINK will attempt to auto-detect available resources upon running regardless of the resources that were requested when the job was submitted.  Therefore, we have to read and parse the requested resources in the cluster config file in order for them to be communicated to PLINK from within the Snakefile.  

### Other notes

It is recommended that snakemake is run as an interactive session on an HPC.  Snakemake will launch the specified number (via the -j flag) of jobs, and then will hang and wait for them to finish.  As jobs finish (and assuming no errors), snakemake will launch additional jobs keeping the total running jobs at whatever -j is set for.  Although Snakemake should not use a lot of memory, it could have long run times, which is generally not advisable on login nodes.  

One attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  To do so, simply resubmit the original snakemake command.

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete

where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

Also, it is often very helpful to do a 'dry-run' of the pipeline in which the different steps and dependencies are printed to screen, but no actual jobs are executed.  This can be helpful to ensure that config entries are correct, etc.  To perform a dry-run, do:

    snakemake -nrp
    
#### Debugging and error reports

Should an error be encountered in a job, snakemake will halt the pipeline and indicate in the terminal that an error has occurred.  The offending job will also be printed in red in the terminal window.  More information on why the job failed can be found in the 'stdout' and 'stderr' files that are output to the **'OandE'** directory and will be labelled with the jobname.

## Pipeline Overview

### Input Data
One or more sets of PLINK files are accepted as input and are specified at the 'data' entry under the 'query' section of the config file.  Along with each dataset, the user can specify the following "keys":  

-   **chrom_key**: tab-delimited text file with 2 columns (no header).  The first column contains the old chromosome names, and the second column contains the new names. 
    - Used for converting to numeric names. e.g chr10 to 10.
-   **allele_key**: tab-delimited text file with 5 columns (no header).  First column is snpID and following columns are: old_allele1 old_allele2 new_allele1 new_allele2.  
    - Used for converting alleles with A/B specification to ACGT.  Oftentimes provided in the dbGaP download.  If alleles are already specified in ACGT format, this field can be set to 'none'
-   **ID_key**: tab-delimited text file with 2 columns (no header).  First column is old SNP ID and second column is new SNP ID.  
    - Used for converting to rsID format.  If SNP IDs are already in rs-format, this field can be set to 'none'
-   **flip_key**: text file with single column containing SNP rsIDs that need to be flipped in order to align strand to the hg19 reference genome. 
    - Used to harmonize strand across datasets to the hg19 reference genome.  Set this field to 'none' if all alleles are already on the same strand as the target reference genome.

Each of these fields are optional and providing 'none' as the entry will disable the steps associated with each key.  However, these fields should only be set to 'none' if you are sure that they are not necessary (e.g. you have already fixed any existing strand issues across datasets).  

Example of input specifications in the config file:

    query:
      "dataset1":
        data: "PATH/TO/PLINK/PREFIX/FOR/DATASET1"
        chrom_key: "PATH/TO/PLINK/PREFIX/FOR/DATASET1"
        allele_key: "PATH/TO/PLINK/PREFIX/FOR/DATASET1"
        ID_key: "PATH/TO/PLINK/PREFIX/FOR/DATASET1"
        flip_key: "PATH/TO/PLINK/PREFIX/FOR/DATASET1"
      "dataset2":
        data: "PATH/TO/PLINK/PREFIX/FOR/DATASET2"
        chrom_key: "none"
        allele_key: "none"
        ID_key: "none"
        flip_key: "none"
        
### Output
The output is a set of PLINK files in the parent directory labelled with the value provided in the 'outname' entry of the config file.  However, if 'merge' is set to 'false' in the config file, this final merge step is skipped, and the final output would be the set of QC'ed plink files within each subdirectory labelled with the dataset names.  Within each of these subdirectories, there will also be a set of VCF files, which are suitable for use in either the Michigan or TOPMed imputation servers.

The other primary output is a PDF report containing a summary of various steps in the pipeline.  It is **highly recommended** that the user carefully review this report to confirm that everything seems in order.  Particular attention should be paid to whether specific steps have resulted in major loss of markers as well as whether there is a positive correlation between allele frequencies in the 1000Genomes dataset and allele frequencies in each of the query datasets.  These scatter plots are provided towards the end of the report, and if a substantial subset of the points exhibit an anti-correlation, this is indicative of a preponderance of strand errors that ought to be corrected (via the 'flip_key') prior to proceeding.     
      
### Dataset harmonization
The first step(s) in the pipeline aims to harmonize the naming of chromosomes, alleles, and variant IDs.  This is accomplished via the 4 keys described above.  While this pipeline generally attempts to simplify the QC process, it is extremely important that the user is acquainted well enough with each individual dataset to ensure that the appropriate keys are specified (or not specified).
   
### Reference allele fixing

In contrast to a VCF, where alleles are specified with respect to a specified reference genome (reference versus alternative alleles), PLINK-formatted files often specify alleles as major/minor alleles based on the frequency in the dataset.  Furthermore, many commonly used arrays will contain a mixture of SNPs genotyped on either the + or - strand.  Lastly, the default behavior of PLINK is to automatically set the minor to A1 and the major allele to A2, which can unintentionally generate inconsistencies in allele specifications across datasets.  

With respect to a reference genome, two possible types of errors can occur:
-   Flipped strand:  The genotype is specified with respect to the opposite strand relative to the reference genome.
-   Swapped allele:  The genotype is specified on the same strand as the reference genome, but the A1 (minor) allele has been set to equal the 'reference' allele when it ought to be set to equal the non-reference/'alternative' allele

To identify these errors, we use the bcftools plugin '+fixref', which requires not only the reference sequence (fasta) file, but also a VCF file containing variant sites that are used to identify mismatching alleles in the query dataset.  Importantly, if the program determines that no strand issues exist and that the reference/alternative alleles have simply been swapped, then program will swap the major/minor alleles to match the reference.  It will not perform any strand flipping, where it converts genotypes to be specified with respect to the nucleotide on the opposite strand.  Although the program will attempt to identify these strand flips, it doesn't make the correction as the authors consider this a risky move that should not be handled in an automated fashion.  Thus, flip-strand mismatches are ultimately removed.  If there are a large number of these, the user should attempt to understand and resolve the source of the issue and rerun this pipeline.

By default, the pipeline will download the following files for the hg19 reference genome:

Reference fasta: 
ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

Reference VCF (1000Genomes):
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz

An indication of whether alleles are now specified correctly is to plot frequency of an allele in the query population against the frequency in the reference population and look for an obviously positive correlation.  Such plots are automatically produced in the PDF report as the final step in the pipeline.

### Basic QC

After alleles have been fixed as described above, a series of basic QC steps are performed on each dataset by the script _'scripts/QC.py'_, with the filtering thresholds specified in the config file (see below).   

    perform_QC: 'true'
    QC:
      vm1: "0.2" # Initial variant missingness filter
      gm: "0.1" # Individual missingness filter
      vm2: "0.05"  # Ultimate call rate for variants after removing low-callrate samples
      maf: "0.01"  # mimimum Minor allele frequency
      hwe: "0.0000001"  # p-value threshold for whether site follows hardy-weinberg
      mbs: "0.0000001"  # p-value treshold for test of whether missingness varies by sex
    
 We first wish to identify and remove individual samples that show high missingess across markers (specified by 'gm').  However, to identify these individuals, we first need to remove variants that imputed poorly across all individuals (specified by 'vm1').  After removing these individuals, we then remove variants with high missingness (specified by 'vm2').  Since poor imputation will result in missing genotypes, this missingness filter indirectly filters for low quality imputation sites.  Variants are also filtered based whether or not they show significant departures from Hardy-Weinberg Equilibrium ('hwe' entry) and whether there is a significant association between missingness and sex ('mbs' entry).  We also remove rare variants based on the 'maf' value.  Lastly, we remove indels, duplicate SNPs, and multi-allelic variants.  
 
 Note that testing for missigness by case/control status is generally recommended as well if the user wishes to proceed straight to SNP-based analyses such as GWAS.  However, if the data is to be used for ancestry inference, it may make more sense to retain these SNPs.  

### Merging inputs (Optional)
If multiple input datasets were provided, an optional final step is to create a single merged dataset consisting of only the sites that overlap (i.e. passed filters) across all component datasets.  This behavior is controlled by the 'merge' entry in the config file.  To enable the merging behavior, set this to:

    merge: 'true'
    
### Imputaton preparation
Another optional, final feature is to create a set of of VCF files (parsed by chromosome) for each of the input datasets.  These VCFs can be used directly as input into either the Michigan Imputation Server or the TOPMed Imputation Server.  The output of the imputation servers can then be used as input into the post-imputation QC pipeline (see README.md in the 'postImpute' directory).