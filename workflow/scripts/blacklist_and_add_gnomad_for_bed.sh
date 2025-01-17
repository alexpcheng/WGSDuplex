#!/bin/bash

input_bed=$1
temp_output_bed=$2
output_bed=$3

gnomad=$4

bedtools subtract -a $input_bed -b resources/hg38/blacklists/encode_blacklist_hg38_v2.bed -header \
| bedtools subtract -a stdin -b resources/hg38/blacklists/ug_lcr.bed -header \
  | bedtools subtract -a stdin -b resources/hg38/blacklists/centromeres.sorted.bed -header \
  | bedtools subtract -a stdin -b resources/hg38/blacklists/simple_repeat_regions.sorted2.bed -header > $temp_output_bed


bedtools intersect -sorted -a $temp_output_bed -b $gnomad -wao > $output_bed


rm $temp_output_bed
