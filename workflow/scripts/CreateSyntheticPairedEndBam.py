#!/usr/bin/env python
# coding: utf-8

import pysam
import pandas as pd
import sys

def rev_comp(seq):
    d = {'A': 'T', 'C':'G', 'G':'C', 'T':'A'}
    seq = seq[::-1]
    rc = ''
    for bp in seq:
        rc+=d[bp]
    return(rc)


def main():
    bam = pysam.AlignmentFile(sys.argv[1], 'rb')
    bam2 = pysam.AlignmentFile(sys.argv[2], 'wb', template = bam)
    i=0

    header = bam.header
    print(header)
    i=0
    for read in bam:
    
        if read.is_reverse:
            R1_flag = 83
            R2_flag = 163
        
        else:
            R1_flag = 99
            R2_flag = 147
    
        if read.is_supplementary:
            R1_flag += 2048
            R2_flag += 2048
        if read.is_secondary:
            R1_flag += 256
            R2_flag += 256

        if read.is_qcfail:
            R1_flag += 512
            R2_flag += 512
        
        RX = read.get_tag('RX')
    
        MC = read.cigarstring
        MQ = read.mapping_quality
    
    
        R1 = pysam.AlignedSegment(header=header)
        R1.query_name = read.query_name
        R1.cigarstring = read.cigarstring
        R1.flag = R1_flag
        R1.query_sequence = read.query_sequence
        R1.query_qualities = read.query_qualities
        R1.reference_id = read.reference_id
        R1.reference_name = read.reference_name
        R1.reference_start = read.reference_start
        R1.mapping_quality = read.mapping_quality
        R1.next_reference_id = read.reference_id
        R1.next_reference_start = R1.reference_start
    
        R1.set_tag('RX', RX)
        R1.set_tag('MC', MC)
        R1.set_tag('MQ', MQ)
    
        R2 = pysam.AlignedSegment(header=header)
        R2.query_name = read.query_name
        R2.cigarstring = read.cigarstring
        R2.flag = R2_flag
        R2.query_sequence = read.query_sequence
        R2.query_qualities = read.query_qualities
        R2.reference_name = read.reference_name
        R2.reference_start = read.reference_start
        R2.mapping_quality = read.mapping_quality
        R2.next_reference_id = read.reference_id
        R2.next_reference_start = R2.reference_start
    
        R2.set_tag('RX', RX)
        R2.set_tag('MC', MC)
        R2.set_tag('MQ', MQ)
    
        bam2.write(R1)
        bam2.write(R2)
        i+=1
    
        if i % 10000000 ==0:
            print(i)
    
    bam.close()
    bam2.close()

if __name__ == "__main__":
    main()
