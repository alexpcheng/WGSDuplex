#!/usr/bin/env python
# coding: utf-8

import pysam
import pandas as pd
import argparse

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
    parser.add_argument('--sscs_bam', help = 'sscs consensus bam file (mapped, sorted and indexed)')
    parser.add_argument('--original_bam', help = 'Original bam, split by chromosome, that we got from UG')
    parser.add_argument('--groupbyumi_bam', help = 'paired-end bam grouped by UMI')
    parser.add_argument('--output_bam', help = 'bam file that only contains reads that duplexed well. Now a SE bam that can go through FlowFeature Mapper')
    
    args = parser.parse_args()
    
    sscs_file = args.sscs_bam
    groupbam_file = args.groupbyumi_bam
    output_file = args.output_bam
    original_file = args.original_bam
    
    return(sscs_file, groupbam_file, output_file, original_file)

def get_sscs_tags(sscs_file):
    i = 0
    sscs_tag_dict = {}
    
    sscs = pysam.AlignmentFile(sscs_file)
    for read in sscs:
        MI = read.get_tag('MI')
        cD = str(read.get_tag('cD'))
        cE = str(read.get_tag('cE'))
        
        cS = str(read.query_sequence)
        zQ = ''.join(chr(x) for x in read.query_qualities) #z because FM sorts alphabetically...
        if read.is_mapped:
            cC = read.cigarstring
            dR = read.reference_name
            dS = read.reference_start+1
        else:
            dR = 'None'
            dS = 'None'
            cC = 'None'
        sscs_tag_dict[MI] = [cD, cE, cS, zQ, dR, dS, cC]
        
        if i % 1000000 == 0:
            print(f'processed {i} sscs reads')
        i+=1
    sscs.close()
    return(sscs_tag_dict)

def get_groupbyumi_counts(groupbam_file):
    i=0
    groupbyumi_count_dict = {}
    readid_to_mi_dict = {}
    
    groupbam = pysam.AlignmentFile(groupbam_file)
    for read in groupbam:
        if read.is_read1:
            MI = read.get_tag('MI')
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


def write_bam(original_file, output_file, sscs_tag_dict, groupbyumi_count_dict, readid_to_mi_dict):
    i=0
    original_bam= pysam.AlignmentFile(original_file)
    new_bam = pysam.AlignmentFile(output_file, 'wb', template = original_bam)
    for read in original_bam:
        qname = read.query_name
    
        if qname in readid_to_mi_dict.keys():
            MI = readid_to_mi_dict[qname]
        else:
            MI = -1
        if MI in groupbyumi_count_dict.keys():
            cU = groupbyumi_count_dict[MI]
        else:
            cU = 0

        if MI in sscs_tag_dict.keys():
            cD, cE, cS, zQ, dR, dS, cC  = sscs_tag_dict[MI]            
            read.set_tag('MI', MI)
            read.set_tag('cU', cU)
            read.set_tag('cD', cD)
            read.set_tag('cE', cE)
            read.set_tag('cS', cS)
            read.set_tag('zQ', zQ)
            read.set_tag('dR', dR)
            read.set_tag('dS', dS)
            read.set_tag('cC', cC)
            read.set_tag('rs', read.reference_start+1)
            new_bam.write(read)

    original_bam.close()
    new_bam.close()
    pysam.index(output_file)
    
def main():
    sscs_file, groupbam_file, output_file, original_file = parse_arguments()
    print("Getting sscs tags")
    sscs_tag_dict = get_sscs_tags(sscs_file)
    print("Getting umi counts")
    groupbyumi_count_dict, readid_to_mi_dict = get_groupbyumi_counts(groupbam_file)
    print("Writting bam")
    write_bam(original_file, output_file, sscs_tag_dict, groupbyumi_count_dict, readid_to_mi_dict)

if __name__ == '__main__':
    main()
