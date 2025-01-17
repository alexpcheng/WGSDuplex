#!/usr/bin/env python
# coding: utf-8

import pysam
import pandas as pd
import argparse
import os
#import Bio
"""
The goal here is to take all the relevant tags from the duplex file that contain 
info that allows us to classify variants as true or errors.

These tags include:
A) cD / aD / bD / cE / MI
These are tags created by fgbio based on the duplex consensus.

B) query sequence and query qualities and reference start position
These values need to be added as tags so FlowFeatureMapper can carry them to 
the final variant VCF INFO column and will be used to get position in read
and BP at that position
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description='Adds relevant tags to a bam file to then be processed by FlowFeatureMapper')
    parser.add_argument('--input_bam')
    parser.add_argument('--output_bam')

    args = parser.parse_args()
    
    input_bam = args.input_bam
    output_bam = args.output_bam
    
    return(input_bam, output_bam)

def get_groupbyumi_counts(groupbam_file):
    i=0
    groupbyumi_count_dict = {}
    readid_to_mi_dict = {}
    
    groupbam = pysam.AlignmentFile(groupbam_file)
    for read in groupbam:
        if read.is_read1:
            MI = read.get_tag('MI').split('/')[0]
            
            readid_to_mi_dict[read.query_name] = MI
            
            if MI in groupbyumi_count_dict.keys():
                groupbyumi_count_dict[MI]+=1
            else:
                groupbyumi_count_dict[MI] = 1
            if i % 1000000 == 0:
                print(f'processed {i} groupedbyumi reads')
            i+=1
    groupbam.close()
    
    return(groupbyumi_count_dict, readid_to_mi_dict)


def write_bam(input_file, output_bam, groupbyumi_count_dict, readid_to_mi_dict):
    i=0
    original_bam= pysam.AlignmentFile(input_file)
    
    new_bam = pysam.AlignmentFile(output_bam+'.tmp.bam', 'wb', template = original_bam)
    
    for read in original_bam:
        if read.is_read1 and not (read.is_supplementary or read.is_secondary):    
            qname = read.query_name

            if qname in readid_to_mi_dict.keys():
                MI = readid_to_mi_dict[qname]
            else:
                MI = -1
            if MI in groupbyumi_count_dict.keys():
                cU = groupbyumi_count_dict[MI]
            else:
                cU = 0

            read.set_tag('MI', MI)
            read.set_tag('cU', cU)
            
            if read.is_forward:
                read.flag = 0
            else:
                read.flag = 16

            new_bam.write(read)

    original_bam.close()
    new_bam.close()
    
    pysam.sort("-o", output_bam, output_bam+".tmp.bam")
    pysam.index(output_bam)
    os.system("rm " + output_bam+".tmp.bam")
    
def main():
    input_bam, output_bam = parse_arguments()

    groupbyumi_count_dict, readid_to_mi_dict = get_groupbyumi_counts(input_bam)
    write_bam(input_bam, output_bam, groupbyumi_count_dict, readid_to_mi_dict)

if __name__ == '__main__':
    main()
