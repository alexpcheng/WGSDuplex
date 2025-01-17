# WGSDuplex


####################################################################################
#   THE FOLLOWING CODE PACKAGE HAS BEEN PREPARED TO PROCESS AND ANALYZE DATA
#   FROM WHOLE GENOME SEQUENCING OF CELL-FREE DNA
#
#
#   THE WORKFLOW CAN BE USED TO GENERATE:
#   1- DUPLEX DENOISED WGS DATA
#   2- SINGLE STRAND CONSENSUS DENOISED WGS DATA
#   3- UMI AGNOSTIC DENOISED WGS DATA
#
#   FOR EACH TYPE OF DENOISING, THE FOLLOWING DATA IS GENERATED
#
#   1- TRINUCLEOTIDE FREQUENCIES OF SOMATIC MUTATIONS
#   2- VARIANTS / GE OF A GIVEN MUTATIONAL PROCESS (EXAMPLE: UV-DERIVED VARIANTS)
####################################################################################

DOWNLOADS

You will need to dowload the following files:
1- human reference genome (hg38; name it "resources/hg38/Homo_sapiens_assembly38.fasta")
2- mouse reference genome and concatenate it to hg38 (name it "resources/hg38_mm39/hg38_mm39.fasta")
3- gnomad pop gen VCF file and subset by chromosome and filter for AF>=0.0001
("resources/hg38/blacklists/gnomad.genomes.v3.1.2.sites.{chrom}.AF0.0001.{bed, vcf.gz}")
4- UG GATK 
Obtain directly from https://www.dropbox.com/scl/fi/r9alrszchfdv2po8wqmvx/gatk-package-ultima_v0.6-SNAPSHOT-local.jar?rlkey=o4aawmvs009ii85ck32wi5mvi&st=k8wetaoo&dl=0
Place it directy in the "software" folder

STEPS TO INSTALL SOFTWARE:

conda env create --name wgsumi --file workflow/envs/wgsumi.yaml
conda env create --name wgsumiR --file workflow/envs/wgsumiR.yaml
conda env create --name AlleleFrequency --file workflow/envs/AlleleFrequency.yaml

conda activate wgsumi
pip install pyfaidx


STEPS TO RUN THE SOFTWARE 
snakemake --use-conda -j {number_of_cores}

EXPECTED OUTPUT
The following files should be created upon a succesful software run:
    'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.duplex_denoised.trinucleotides.txt'
    'results/WgsDuplex/SingleStrand/{sample_id}/{sample_id}.singlestrand_denoised.trinucleotides.txt'
    'results/WgsDuplex/UmiAgnostic/{sample_id}/{sample_id}.umiagno_denoised.trinucleotides.txt'


VERSIONING
Software has been test on the following OS:
    Operating System: CentOS Linux 7 (Core)
    CPE OS Name: cpe:/o:centos:centos:7
    Kernel: Linux 3.10.0-1160.118.1.el7.x86_64
    Architecture: x86-64

Using a 64GB, 16 core machine, the pipeline took approximately X hours to complete on the provided dataset


