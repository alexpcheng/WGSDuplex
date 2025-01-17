#!/usr/bin/env python
# coding: utf-8


import pysam
import pandas as pd
import argparse
import re
import sys


"""
We are going to take the FeatureMap file and do the following:
1- Add duplex base and position in read (PIR) and position from end of read (PIR2)
2- Add lofreq statistics (depth, af, sb)
3- How many reads passed through the FeatureMap filters
4- Quality of duplex (how many Ns)
"""

def parse_arguments():
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    
    return(in_file, out_file)    

def pass_line(line):
    chrom, pos, ref, alt, altID, RX, \
    X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT,\
    X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE, rq, rs, \
    PIR, PIR2, \
    DP, AF, SB, ref_top, ref_bot, alt_top, alt_bot, gnomad_chrom, gnomad_pos, gnomad_id, gnomad_ref, gnomad_alt = line.strip().split('\t')
    
    p = True
    
    PIR = int(PIR)
    PIR2= int(PIR2)
    
    if not (X_FLAGS == "0" or X_FLAGS == "16"):
        p = False
        
    if X_FLAGS == "0" and alt_top !="None":
        if int(alt_top) < 1:
            p = False
    if X_FLAGS == "0" and alt_top =="None": 
        p = False
        
    if X_FLAGS == "16" and alt_bot =="None": 
        p = False
            
    if X_FLAGS == "16" and alt_bot !="None":
        if int(alt_bot) < 1:
            p = False
    
    if PIR < 0 or PIR2 <0:
        p = False
    if PIR <=10 or PIR2 <=10:
        p = False
        
    ## Extra filtering to do our best...
    if int(X_LENGTH) >= 200:
        p = False
    if int(X_MAPQ) < 60:
        p = False
    if int(X_FC1) >=10:
        p = False
    if int(X_EDIST) >4:
        p = False
    
    newline = '\t'.join([chrom, pos, ref, alt, altID, RX,\
    X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT,\
    X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE, rq, rs, \
    str(PIR), str(PIR2), \
    DP, AF, SB, ref_top, ref_bot, alt_top, alt_bot]) + '\n'
    
    return(newline, p)
    
def filter_featuremap(in_file, out_file):
    with open(in_file) as f, open(out_file, 'w') as w:
        next(f)
        w.write('\t'.join(['chrom', 'pos', 'ref', 'alt', 'altID', 'RX',\
                          'X_CIGAR', 'X_EDIST', 'X_FC1', 'X_FC2', 'X_FILTERED_COUNT',\
                          'X_FLAGS', 'X_LENGTH', 'X_MAPQ', 'X_READ_COUNT', 'X_RN', 'X_SCORE', 'rq', 'rs', \
                          'PIR', 'PIR2', 'DP', 'AF', 'SB', 'ref_top', 'ref_bot', 'alt_top', 'alt_bot']) + '\n')
        for line in f:
            newline, p = pass_line(line)
            if p:
                w.write(newline)

def main():
    in_file, out_file = parse_arguments()
    
    filter_featuremap(in_file, out_file)


if __name__ == '__main__':
    main()

