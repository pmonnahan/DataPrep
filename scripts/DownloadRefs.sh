
# Download all the Reference Data to Reformat the files
# ----------------------------------------------------------------------------

if [ ! -f ./RefAnnotationData/human_g1k_v37.fasta.gz ]; then

echo
echo Downloading Reference Data and index files from 1K Genomes and NCBI
echo ----------------------------------------------

echo Downloading: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
echo
#wget --directory-prefix=./RefAnnotationData/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#wget --directory-prefix=./RefAnnotationData/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

fi

# Download the annotation files (make sure the the build version is correct) to flip/fix the alleles
if [ ! -f ./RefAnnotationData/All_20170710.vcf.gz ]; then
echo
echo Downloading ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
echo
wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi


fi

# For tracking completion via snakemake
echo "downloaded" > ./RefAnnotationData/downloaded.txt