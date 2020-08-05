# Post-imputation QC pipeline

The purpose of this pipeline is to convert(/combine) one or more imputed datasets created by the [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/) into a single set (separated by chromosome) of PLINK-formatted datasets.  In brief, one or more sets of gzipped, _human_ VCF files are first run through [CrossMap](http://crossmap.sourceforge.net/), which converts coordinates from the GRCh38 reference genome to GRCh19.  Then, these files are converted to PLINK format and variants are filtered for missingness, duplicates, and indels (are removed).  For each chromosome, we then merge these resulting files across datasets.  Only variants that have been retained across all datasets are included in this merged dataset.  Rare alleles are then filtered from this merged dataset.  

![Pipeline DAG](https://github.com/pmonnahan/DataPrep/blob/master/postImpute/Pipeline_DAG.png)



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
   - [Coordinate Conversion](#coordinate-conversion)
   - [PLINK Conversion and initial QC](#plink-conversion-and-initial-qc)
   - [Merging Datasets](#merging-datasets)
   - [Optional extraction of specific SNPS](#optional-extraction-of-specific-snps)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Requirements

### Snakemake
The pipeline is coordinated and run on an HPC (or locally) using _Snakemake_.  On UMN's HPC, snakemake can be installed by:

    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge -c bioconda snakemake python=3.6

The 'module load' command will likely need to be run each time prior to use of Snakemake.

Alternatively, you can try installing snakemake via _pip_:

    pip3 install --user snakemake pyaml

### Singularity

The installation of the individual programs used throughout this pipeline can be completely avoid by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image [here](https://github.com/pmonnahan/AncInf/blob/master/singularity/Singularity_defs.def).  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

However, in order to utilize the singularity image, _Singularity_ must be installed on the HPC.  Currently, the pipeline assumes that _Singularity_ will be available as a module and can be loaded into the environment via the command specified in the config.yml file, where it says 'singularity_module'.  The default setting will work for MSI at UMN.

Singularity settings in config.yml

    singularity:
      use_singularity: 'true'
      image: '/home/pmonnaha/pmonnaha/singularity/AncestryInference.sif


## Running the workflow

Clone the parent repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/DataPrep.git postImputeQC
    
Then, do

    cd postImputeQC
    rm -r workflow
    mv postImpute/workflow .
    
    
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
The pipeline accepts as input one or more directories containing the gzipped VCFs from the TOPMed imputation server.  Each of these directories should be listed under the 'query' section of the config file as demonstrated below.

    query:
      "dataset1": "/PATH/TO/DATASET1/DIRECTORY"
      "dataset2": "/PATH/TO/DATASET2/DIRECTORY"
      "dataset3": "/PATH/TO/DATASET3/DIRECTORY"

Additionally, these input directories should contain two additional files that are provided by the imputation server called 'chunks-excluded.txt' and 'snps-excluded.txt'.  If you submitted separate jobs for each chromosome, then these files will need to be concatenated and provided as a single pair of files.

Phenotypes of the samples must be specified in one of two ways.  When more than one input dataset is provided, if the individual datasets are composed entirely of cases OR controls, then the names of the datasets as listed under 'query' in the config.yml can simply be specified under the 'case_datasets' and 'control_dataset' fields in the config file, respectively (see below).  If multiple case or multiple control datasets exist, these can be listed as a comma-separated (no spaces) string in these fields.  The other way to specify phenotypes is to provide a tab-delimited text file where the first column contains the sample IDs (as they appear in the imputed VCF file) and the second column contains the phenotype. The path to this file can be provided in the field labelled 'phenotype_file' under the 'phenotype_data' field in the config.yml file.

Sex of the samples must also be specified in a tab-delimited text file where the first column is sample ID and the second column is the sex specification according to PLINK.  The path to this file can be provided in the field labelled 'sex_file' under the 'phenotype_data' field in the config.yml file.

    phenotype_data: 
      pheno_file: "none"
      case_datasets: "dataset1,dataset2"
      control_datasets: "dataset3"
      sex_file: "/path/to/sex/file"
      
### Output
The final output of this pipeline is a set of PLINK files which have been named according to the entry provided in the 'outname' of the config file, e.g.:

    outname: "Example"

Additionally, a PDF report bearing this prefix is created, which summarizes much of the relevant information from the different steps of the pipeline.  
 
### Coordinate Conversion

The first step in the pipeline is to convert coordinates from GRCh38 to GRCh19 using the program [CrossMap](http://crossmap.sourceforge.net/).  The default behavior, as indicated in the config settings listed below, is to download the reference fasta and 'chain' files (key linking coordinates across chromosomes).  If you wish to download a different set of files, the links in each entry can be modified.  Alternatively, if you wish to use files that you have already procured, set the 'download_refs' entry to 'false', and provide the paths in the 'chain' and 'fasta' entries.  Note that these files must be filtered to exclude any non-canonical chromosome names.

    CrossMap:
      download_refs: 'true'
      chain: "ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz"  
      fasta: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
      threads: "6"
      exec: "/usr/local/bin/miniconda/bin/CrossMap.py"
      
### PLINK Conversion and initial QC
Following coordinate conversion, the VCFs for each chromosome are converted to plink format.  During this conversion, poorly imputed genotypes are set to missing.  That is, if the probability of the most probable genotype falls below the threshold specified in the config file, then the genotype for this sample is set to missing.  The relevant entry in the config file is:
      
    min_gp: "0.85" 
    
Then, the variant IDs are updated to reflect the new coordinates, and the sex of the samples is incorporated into the PLINK files, using the 'sex_file' provided in the config file.  
    
A series of filtering steps are then performed by the script _'scripts/QC.py'_, with the filtering thresholds specified in the config file (see below).  

    perform_QC: 'true'
    QC:
      vm1: "0.2" # Initial variant missingness filter
      gm: "0.1" # Individual missingness filter
      vm2: "0.05"  # Ultimate call rate for variants after removing low-callrate samples
      maf: "0.01"  # mimimum Minor allele frequency
      hwe: "0.0000001"  # p-value threshold for whether site follows hardy-weinberg
      mbs: "0.0000001"  # p-value treshold for test of whether missingness varies by sex
    
 We first wish to identify and remove individual samples that show high missingess across markers (specified by 'gm').  However, to identify these individuals, we first need to remove variants that imputed poorly across all individuals (specified by 'vm1').  After removing these individuals, we then remove variants with high missingness (specified by 'vm2').  Since poor imputation will result in missing genotypes, this missingness filter indirectly filters for low quality imputation sites.  Variants are also filtered based whether or not they show significant departures from Hardy-Weinberg Equilibrium ('hwe' entry) and whether there is a significant association between missingness and sex ('mbs' entry).  Note that the minor allele frequency filter ('maf') is reserved for a later step following the (optional) merging of input datasets.  
  
  Lastly, we remove indels, duplicate SNPs, and multi-allelic variants.  Note that testing for missigness by case/control status is generally recommended as well, but is implemented early on in a downstream pipeline for association/admixture mapping (i.e [admixMap](https://github.com/pmonnahan/admixMap)).  
 
 At this point, phenotypes are incorporated into the PLINK files, and the reference/alternative alleles are re-specified as in the raw plink files.  While the correct allele specification should have been retained throughout, we include this step here just in case any unexpected allele swapping occurred during QC.  
 
### Merging Datasets
If multiple imputed datasets were provided as input, the next step would be, for each chromosome, to merge the genotypes across datasets.  Importantly, only variants that are still present in all datasets (i.e. have not been filtered in any single dataset) will be retained.  This way, if a variant imputed poorly in one dataset for whatever reason, it would be removed entirely from the merged dataset.  

The merged datasets are then filtered for minor allele frequency, whose value is specified in the 'QC' section of the config file discussed above.  This filter is reserved for this step, such that alleles that are rare in one dataset, but common in another, are not lost prior to merging.  
  
### Optional extraction of specific SNPS
A user may wish to ensure that certain SNPs are retained throughout this process, regardless of whether they pass all filters.  For example, a prior GWAS reported effects at specific SNPs that the user wants to include in a subsequent analysis, but one or more of these SNPs imputed poorly.  A tab-delimited file containing the chromosome and position (one SNP per line) can be provided at the 'snp_list' entry.  

    snp_list: "accessory/ALL_snpLoc.txt"
    
After the pipeline has completed, these specific SNPs can be extracted by running:

    snakemake <outname>.snps.bed --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 1
    
, where <outname> should be replaced by the value used in the 'outname' entry of the config file.  
