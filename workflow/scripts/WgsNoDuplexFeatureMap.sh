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
#chrom=$3

ref=/gpfs/commons/projects/high_depth_liquid_biopsy/parameter_data/Homo_sapiens_assembly38.fasta
intervals=/gpfs/commons/projects/high_depth_liquid_biopsy/parameter_data/wgs_calling_regions.hg38.interval_list
encode=/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/PROBLEMATIC_REGIONS/encode_blacklist_hg38_v2.bed

ug_lcr=/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/PROBLEMATIC_REGION_ANALYSIS/ug_lcr.bed
centro=/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/PROBLEMATIC_REGION_ANALYSIS/centromeres.sorted.bed
simple=/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/PROBLEMATIC_REGION_ANALYSIS/simple_repeat_regions.sorted2.bed
gnomad=/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/PROBLEMATIC_REGION_ANALYSIS/gnomad.genomes.v3.1.2.sites.autosomalXY.AF0.0001.bed

java -Xmx80G -jar /gpfs/commons/projects/high_depth_liquid_biopsy/gatk-private_ultima_v0.6/gatk-private/build/libs/gatk-package-ultima_v0.6-SNAPSHOT-local.jar FlowFeatureMapper \
	-I $bam \
	-O $output \
	-R $ref \
	-XL $encode \
	-XL $ug_lcr \
	-XL $centro \
	-XL $simple \
	-XL $gnomad \
	--intervals $intervals \
	--snv-identical-bases 5 \
	--snv-identical-bases-after 5 \
	--min-score 3 \
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
	--copy-attr rs
