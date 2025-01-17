#!/usr/bin/env python
# coding: utf-8


import pysam
import pandas as pd
import argparse
import re
import glob
from natsort import natsorted
from pyfaidx import Fasta

"""
We are going to take the FeatureMap file and do the following:
1- Do our best
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description='Add a couple features you know how it is')
    parser.add_argument('--sample_id')
    parser.add_argument('--path')
    parser.add_argument('--output')
    
    args = parser.parse_args()
    
    sample_id = args.sample_id
    path = args.path
    output = args.output
   
    
    return(sample_id, path, output)


def cigarstring_to_tuple(cigarstring):
    cigarlist = list(filter(None, re.split(r'(\d+)', cigarstring)))
    cigar = []
    lengths = cigarlist[::2] 
    operations = cigarlist[1::2]

    
    pysam_dict = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8, 'B':9}
    for op, le in zip(operations, lengths):
        op2 = pysam_dict[op]
        cigar.append((op2, int(le)))
    return(cigar)
    
    
    
def get_PIR(cC, position, dS):
    """ Obtains the correct position in read of a given variant position. Accounts for soft-clipping and indels.
    read: read from bam
    position: coordinate of variant
    """
    start = dS-1 #to confirm to pysam
    cigar = cigarstring_to_tuple(cC)
    if cigar[0][0] == 4:
        start -= cigar[0][1]
    pir_raw = int(position) - start - 1

    if len(cigar) == 1 and cigar[0][0] == 0:
        return pir_raw

    pir = pir_raw
    cigar_counter = 0
    for operator, length in cigar:
        if pir >= cigar_counter + length:
            if operator == 1:
                pir += length
            elif operator == 2:
                pir -= length
        elif pir < cigar_counter:
            return pir
        else:
            if operator == 1:
                pir += length
            elif operator in [2, 4]:
                return -1
        if operator != 2:
            cigar_counter += length
    return pir

def parse_lines(line, hg38):
    chrom, pos, ref, alt, altID, MI, MI_pos, RX, X_CIGAR, X_EDIST, \
    X_FC1, X_FC2, X_FILTERED_COUNT, X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE, rq, \
    rs, aD, aE, bD, bE, cD, cE, cU, cS, dR, \
    dS, cC, FM_passed, num_N, duplex_length, duplex_PIR, duplex_PIR2, duplex_bp, DP, AF, \
    SB, ref_top, ref_bot, alt_top, alt_bot, A_top, A_bot, T_top, T_bot, C_top, \
    C_bot, G_top, G_bot, N_top, N_bot, duplex_bp2 = line.strip().split('\t')
    pir = get_PIR(X_CIGAR, int(pos), int(rs))
    pir2 = int(X_LENGTH) - pir
    
    triN = str(hg38[chrom][int(pos)-2]).upper() + ref + str(hg38[chrom][int(pos)]).upper()
#     print(str(hg38[chrom][int(pos)-1]).upper())
#     print(ref)
#     print(triN)
#     print("==========")
    
    line2 = line.strip()
    line2 = line2 + '\t' + str(pir) + '\t' + str(pir2) + '\t' + triN + '\n'
    return(line2)

def main():
    sample_id, path, output = parse_arguments()
    
    hg38 = Fasta('resources/hg38/Homo_sapiens_assembly38.fasta')
    
    head = True
    w = open(output, 'w')
    
    files = glob.glob(f'{path}/*/*annotatedfeaturemap.mlduplex.tsv')
    files = natsorted(files)
    for file in files:
        with open(file) as f:
            for line in f:
                if head:
                    print(line)
                    w.write(line.strip() + '\t' + 'pir' + '\t' + 'pir2' + '\t' + 'triN' + '\n')
                    head = False
                else:
                    if "chrom" not in line:
                        w.write(parse_lines(line, hg38))
    w.close()
        


if __name__ == '__main__':
    main()

