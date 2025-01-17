#!/bin/bash
#SBATCH --job-name=FeatureMap
#SBATCH --partition=pe2
#SBATCH --mail-user=alcheng@nygenome.org
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=80G
#SBATCH --time=167:00:00
#SBATCH --output=stdout_%j.log
#SBATCH --cpus-per-task=1

module load java

bam=$1
output=$2
chrom=$3

ref=resources/hg38/Homo_sapiens_assembly38.fasta
intervals=resources/hg38/wgs_calling_regions.hg38.interval_list
encode=resources/hg38/blacklists/encode_blacklist_hg38_v2.bed

ug_lcr=resources/hg38/blacklists/ug_lcr.bed
centro=resources/hg38/blacklists/centromeres.sorted.bed
simple=resources/hg38/blacklists/simple_repeat_regions.sorted2.bed
gnomad=resources/hg38/blacklists/gnomad.genomes.v3.1.2.sites.$chrom.AF0.0001.bed

java -Xmx30G -jar software/gatk-package-ultima_v0.6-SNAPSHOT-local.jar FlowFeatureMapper \
	-I $bam \
	-O $output \
	--include-dup-reads \
	-R $ref \
	-XL $encode \
	-XL $ug_lcr \
	-XL $centro \
	-XL $simple \
	-XL $gnomad \
	--intervals $intervals \
	--snv-identical-bases 5 \
	--snv-identical-bases-after 5 \
	--min-score 0 \
	--limit-score 10 \
	--read-filter MappingQualityReadFilter \
	--minimum-mapping-quality 60 \
	--smith-waterman FASTEST_AVAILABLE \
	--likelihood-calculation-engine FlowBased \
	-mbq 0 \
	--kmer-size 10 \
	--copy-attr tr \
	--copy-attr tt \
	--copy-attr rq \
	--copy-attr RX \
	--copy-attr MI \
	--copy-attr cU \
	--copy-attr cD \
	--copy-attr cE \
	--copy-attr cS \
	--copy-attr dR \
	--copy-attr dS \
	--copy-attr cC \
	--copy-attr rs
