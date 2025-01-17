#!/usr/bin/env python
# coding: utf-8


import pysam
import pandas as pd
import argparse
import re


"""
We are going to take the FeatureMap file and do the following:
1- Add duplex base and position in read (PIR) and position from end of read (PIR2)
2- Add lofreq statistics (depth, af, sb)
3- How many reads passed through the FeatureMap filters
4- Quality of duplex (how many Ns)
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description='Add a couple features you know how it is')
    parser.add_argument('--FeatureMap', help = 'FeatureMap file with the necessary tags')
    parser.add_argument('--lofreq_vcf', help = 'lofreq file')
    parser.add_argument('--output')
    parser.add_argument('--duplex_only_output')
    
    args = parser.parse_args()
    
    featuremap_file = args.FeatureMap
    output_file = args.output
    dupout_file = args.duplex_only_output
    lofreq_file = args.lofreq_vcf
    
    return(featuremap_file, output_file, lofreq_file, dupout_file)


def split_featuremap_line(line):
    #print(line)
    chrom, pos, ID, ref, alt, a, b, INFO = line.strip().split('\t')

    x = [i.split('=')[1] for i in INFO.split(';')[0:26]]
    #Need to split at 22 because query qualities can contain ';' and '=' characters
    #zQ= INFO.split('zQ=')[1] # duplex qualities problem for another day
    
    MI = x[0]                 # MI = Molecular ID
    RX = x[1]                 # RX = UMI tag
    X_CIGAR = x[2]            # X_CIGAR = cigar string of the uncollapsed read
    X_EDIST = x[3]            # X_EDIST = edit distance of the uncollapsed read
    X_FC1 = x[4]              # X_FC1 = number of features (SNPs) on the uncollapsed read (pre-filtering)
    X_FC2 = x[5]              # X_FC2 = number of features (SNPs) on the uncollapsed read (post-filtering)
    X_FILTERED_COUNT = x[6]   # X_FILTERED_COUNT = coverage at that position (post-filtering)
    X_FLAGS = x[7]            # X_FLAGS = mapping flag of uncollapsed read
    X_LENGTH = x[8]           # X_LENGTH = molecule length
    X_MAPQ = x[9]             # X_MAPQ = mapping quality
    X_READ_COUNT = x[10]      # X_READCOUNT = coverage at that position (pre-filtering)
    X_RN = x[11]              # X_RN = read id
    X_SCORE = x[12]           # X_SCORE = flow score
    aD = x[13]                # aD = depth of top strands for duplex
    aE = x[14]                # aE = error rate of top strand consensus sequence
    bD = x[15]                # bD = depth of bot strands for duplex
    bE = x[16]                # bE = error rate of bot strand consensus sequence
    cC = x[17]
    cD = x[18]                # cD = depth of the duplex
    cE = x[19]                # cE = error rate of the duplex
    cS = x[20]                # cS = consensus sequence of the duplex
    cU = x[21]                # cU = total number of reads with the duplex UMI (this might not be = cD if some reads had inconsistent cigars)
    dR = x[22]                # dR = chromosme of the duplex
    dS = x[23]                # dS = read start of the duplex
    rq = x[24]
    rs = x[25]
    return(chrom, pos, ref, alt, MI, RX, \
               X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT, X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE, \
               aD, aE, bD, bE, cC, cD, cE, cU, cS, dR, dS, rq, rs)


# In[ ]:





# In[147]:


def count_MIs_that_pass_FM(featuremap_file):
    FM_passed_MI = {}
    with open(featuremap_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            chrom, pos, ref, alt, MI, RX,\
            X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT, X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE,\
            aD, aE, bD, bE, cC, cD, cE, cU, cS, dR, dS, rq, rs = split_featuremap_line(line)

            MI_pos = MI + '|' + chrom + '|' + pos
        
            if MI_pos not in FM_passed_MI.keys():
                FM_passed_MI[MI_pos] = 1
            else:
                FM_passed_MI[MI_pos] +=1
    return(FM_passed_MI)


# In[168]:


def get_lofreq_features(lofreq_file):
    lofreq_feats = {}
    with open(lofreq_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, pos, ID, ref, alt, x, y, INFO = line.strip().split('\t')
            DP = INFO.split(';')[0].split('=')[1]
            AF = INFO.split(';')[1].split('=')[1]
            SB = INFO.split(';')[2].split('=')[1]
            DP4 = INFO.split(';')[3].split('=')[1]
            ref_top = DP4.split(',')[0]
            ref_bot = DP4.split(',')[1]
            alt_top = DP4.split(',')[2]
            alt_bot = DP4.split(',')[3]
            altID = f'{chrom}|{pos}|{ref}|{alt}'
            lofreq_feats[altID] = [DP, AF, SB, ref_top, ref_bot, alt_top, alt_bot]
    return(lofreq_feats)


# In[188]:


def annotate_FM(featuremap_file, output_file, FM_passed_MI, lofreq_feats, dupout_file):
    with open(featuremap_file) as f, open(output_file, 'w') as w, open(dupout_file, 'w') as w2:
        w.write('\t'.join(['chrom', 'pos', 'ref', 'alt', 'altID', 'MI', 'MI_pos', 'RX',\
                           'X_CIGAR', 'X_EDIST', 'X_FC1', 'X_FC2', 'X_FILTERED_COUNT',\
                           'X_FLAGS', 'X_LENGTH', 'X_MAPQ', 'X_READ_COUNT', 'X_RN', 'X_SCORE', 'rq', 'rs', \
                           'aD', 'aE', 'bD', 'bE', 'cD', 'cE', 'cU', 'cS', 'dR', 'dS', 'cC', 'FM_passed',\
                           'num_N', 'duplex_length', 'duplex_PIR', 'duplex_PIR2', 'duplex_bp', \
                           'DP', 'AF', 'SB', 'ref_top', 'ref_bot', 'alt_top', 'alt_bot'])+'\n')
        w2.write('\t'.join(['chrom', 'pos', 'ref', 'alt', 'altID', 'MI', 'MI_pos', 'RX',\
                           'X_CIGAR', 'X_EDIST', 'X_FC1', 'X_FC2', 'X_FILTERED_COUNT',\
                           'X_FLAGS', 'X_LENGTH', 'X_MAPQ', 'X_READ_COUNT', 'X_RN', 'X_SCORE', 'rq', 'rs', \
                           'aD', 'aE', 'bD', 'bE', 'cD', 'cE', 'cU', 'cS', 'dR', 'dS', 'cC', 'FM_passed',\
                           'num_N', 'duplex_length', 'duplex_PIR', 'duplex_PIR2', 'duplex_bp', \
                           'DP', 'AF', 'SB', 'ref_top', 'ref_bot', 'alt_top', 'alt_bot'])+'\n')
        for line in f:
            if line.startswith('#'):
                continue
                
            chrom, pos, ref, alt, MI, RX,\
            X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT, X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE,\
            aD, aE, bD, bE, cC, cD, cE, cU, cS, dR, dS, rq, rs= split_featuremap_line(line)
            
            MI_pos = MI+ '|' + chrom + '|' + pos
            
            FM_passed = FM_passed_MI[MI_pos]
            
            altID = f'{chrom}|{pos}|{ref}|{alt}'
            
            if altID in lofreq_feats.keys():
                DP, AF, SB, ref_top, ref_bot, alt_top, alt_bot = lofreq_feats[altID]
            else:
                DP = 'None'
                AF = 'None'
                SB = 'None'
                ref_top = 'None'
                ref_bot = 'None'
                alt_top = 'None'
                alt_bot = 'None'
            
            if cS != 'None':
                num_N = cS.count('N')
                duplex_length = len(cS)
                if dR == chrom and (int(pos) >= int(dS) and int(pos) < int(dS)+duplex_length):
                    duplex_PIR = get_PIR(cS, cC, int(pos), int(dS))
                    if duplex_PIR < duplex_length:
                        duplex_PIR2 = duplex_length - duplex_PIR
                        duplex_bp = cS[duplex_PIR]
                    else:
                        duplex_PIR = 'None'
                        duplex_PIR2 = 'None'
                        duplex_bp = 'None'
                else:
                    duplex_PIR = -1
                    duplex_PIR2 = -1
                    duplex_bp = 'None'
            else:
                num_N = -1
                duplex_length = -1
                duplex_PIR = -1
                duplex_PIR2 = -1
                duplex_bp = 'None'
            
            line_list = [chrom, pos, ref, alt, altID, MI, MI_pos, RX,\
                         X_CIGAR, X_EDIST, X_FC1, X_FC2, X_FILTERED_COUNT,\
                         X_FLAGS, X_LENGTH, X_MAPQ, X_READ_COUNT, X_RN, X_SCORE, rq, rs, \
                         aD, aE, bD, bE, cD, cE, cU, cS, dR, dS, cC, FM_passed,\
                         num_N, duplex_length, duplex_PIR, duplex_PIR2, duplex_bp,\
                         DP, AF, SB, ref_top, ref_bot, alt_top, alt_bot]
            line_list2 = [str(x) for x in line_list]
            w.write('\t'.join(line_list2)+'\n')
            if cS != 'None':
                w2.write('\t'.join(line_list2)+'\n')
        

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
    
    
    
def get_PIR(cS, cC, position, dS):
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


def main():
    featuremap_file, output_file, lofreq_file, dupout_file = parse_arguments()
    FM_passed_MI = count_MIs_that_pass_FM(featuremap_file)


    lofreq_feats = get_lofreq_features(lofreq_file)


    annotate_FM(featuremap_file, output_file, FM_passed_MI, lofreq_feats, dupout_file)


if __name__ == '__main__':
    main()

