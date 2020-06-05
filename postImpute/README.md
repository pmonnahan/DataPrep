#Post-imputation QC pipeline

The purpose of this pipeline is to convert(/combine) one or more imputed datasets created by the [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/) into a single set (separated by chromosome) of PLINK-formatted datasets.  In brief, one or more sets of gzipped, _human_ VCF files are first run through [CrossMap](http://crossmap.sourceforge.net/), which converts coordinates from the GRCh38 reference genome to GRCh19.  Then, these files are converted to PLINK format and variants are filtered for missingness, duplicates, and indels (are removed).  For each chromosome, we then merge these resulting files across datasets.  Only variants that have been retained across all datasets are included in this merged dataset.  Rare alleles are then filtered from this merged dataset.  

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


## Data preparation

It is assumed that all of the gzipped VCF files are located in a single input directory (specified in the config.yml file; see below).  Additionally, this input directory should contain two additional files that are provided by the imputation server called 'chunks-excluded.txt' and 'snps-excluded.txt'.  If you submitted separate jobs for each chromosome, then these files will need to be concatenated and provided as a single pair of files.

Phenotypes of the samples can be specified in one of two ways.  When more than one input dataset is provided, if the individual datasets are composed entirely of cases OR controls, then the names of the datasets as listed under 'query' in the config.yml can simply be specified under the 'case_datasets' and 'control_dataset' fields in the config file, respectively.  If multiple case or multiple control datasets exist, these can be listed as a comma-separated (no spaces) string in these fields.  The other way to specify phenotypes is to provide a tab-delimited text file where the first column contains the sample IDs (as they appear in the imputed VCF file) and the second column contains the phenotype. The path to this file can be provided in the field labelled 'phenotype_file' under the 'phenotype_data' field in the config.yml file.

Sex of the samples can also be specified in a tab-delimited text file where the first column is sample ID and the second column is the sex specification according to PLINK.  The path to this file can be provided in the field labelled 'sex_file' under the 'phenotype_data' field in the config.yml file.

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