#Pre-imputation QC pipeline

The purpose of this pipeline is to perform the following for a set of input PLINK datasets:

-   basic QC (genotype/variant missingness, HWE, and minor allele frequency)
-   harmonize allele specifications with the GRCh37 reference genome
-   produce a set of VCF files (separated by chromosome) for imputation
-   merge filtered datasets into a single dataset consisting only of overlapping sites

## Data preparation

Following the 'data' field (containing the path info) of each query dataset specified in the config.yml file, the following four fields must be specified.

-   **chrom_key**: tab-delimited text file with 2 columns (no header).  The first column contains the old chromosome names, and the second column contains the new names. Used for converting to numeric names. e.g chr10 to 10.
-   **allele_key**: tab-delimited text file with 5 columns (no header).  First column is snpID and following columns are: old_allele1 old_allele2 new_allele1 new_allele2.  Used for converting alleles with A/B specification to ACGT.  Oftentimes provided in the dbGaP download.  If alleles are already specified in ACGT format, this field can be set to 'none'
-   **ID_key**: tab-delimited text file with 2 columns (no header).  First column is old SNP ID and second column is new SNP ID.  Used for converting to rsID format.  If SNP IDs are already in rs-format, this field can be set to 'none'
-   **flip_key**: text file with single column containing SNP rsIDs that need to be flipped. Set this field to 'none' if all alleles are already on the same strand as the target reference genome.

For an example of how all components for a query dataset could be specified:

    query:
      "cog1":
        data: "/home/pmonnaha/shared/data/admiral/raw/phg000408.v1.NIGMS_ALL_Relapse_GWAS.genotype-calls-
    matrixfmt.c1.HMB/with_fam/cog9904_9905_snp6_with_fam"
        chrom_key: "/home/pmonnaha/pmonnaha/misc/PlinkChrRename.txt"
        allele_key: "/home/pmonnaha/shared/data/admiral/raw/phg000408.v1.NIGMS_ALL_Relapse_GWAS.marker-in
    fo.c.MULTI/recodePlinkAlleles.txt"
        ID_key: "/home/pmonnaha/shared/data/admiral/raw/phg000408.v1.NIGMS_ALL_Relapse_GWAS.marker-info.c
    .MULTI/recodePlinkMarkers_noDups.txt"
        flip_key: "/home/pmonnaha/shared/data/admiral/raw/phg000408.v1.NIGMS_ALL_Relapse_GWAS.marker-info
    .c.MULTI/opposite_strands_rsIDs.txt"

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

The **Snakefile** is the primary workhouse of _snakemake_, which specifies the dependencies of various parts of the pipeline and coordinates their submission as jobs to the MSI cluster.  No modifications to the **Snakefile** are necessary.  However, in order for the **Snakefile** to locate all of the necessary input and correctly submit jobs to the cluster, **both** the config.yaml and cluster.yaml need to be modified.  Open these files and change the entries that are indicated with 'MODIFY'.  

Once these files have been modified, the entire pipeline can be run from within the cloned folder via:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32

where -j specifies the number of jobs that can be submitted at once.  

Additional features regarding the HPC, such as queue and requested resources, are specified in the workflow/cluster.yaml file.

The attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete

where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.