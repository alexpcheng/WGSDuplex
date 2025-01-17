import os
import pysam
import argparse
import time
import pandas as pd
import numpy as np

"""
FILTERING CRITERIA

X_LENGTH<=200
MAPQ==60
X_SCORE>=3
EDIT DISTANCE <=3
"""
def parse_arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_bam")
    parser.add_argument("--output_bam")
    parser.add_argument("--output_featuremap")
    parser.add_argument("--tumor_vcf", default = None)
    parser.add_argument("--output_tf")
    parser.add_argument("--sample_id")
    
    args = parser.parse_args()
    
    input_bam = args.input_bam
    output_bam = args.output_bam
    output_featuremap = args.output_featuremap
    tumor_vcf = args.tumor_vcf
    output_tf = args.output_tf
    sample_id = args.sample_id
    
    arguments = [input_bam, output_bam, output_featuremap, tumor_vcf, output_tf, sample_id]
    
    return(arguments)

def intersect_tumor_bam(input_bam, output_bam, tumor_vcf):
    
    bam = pysam.AlignmentFile(input_bam, 'rb')
    tmp_bam = output_bam.replace(".bam", ".tmp.bam")
    outbam = pysam.AlignmentFile(tmp_bam, 'wb', template = bam)
    read_ids = {}
    i=1
    wrote=0
    t=time.time()
    with open(tumor_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                chrom, pos = line.strip().split('\t')[0:2]
                pos = int(pos)
                for read in bam.fetch(chrom, pos-5, pos+5):
                    if read.query_name not in read_ids.keys():
                        X_LENGTH=len(read.query_sequence)
                        MAPQ = read.mapping_quality
                        edist = read.get_tag('NM')
                        if read.is_secondary or read.is_supplementary or read.is_duplicate or X_LENGTH>200 or MAPQ<60 or edist >=4:
                            continue
                        else:
                            outbam.write(read)
                            wrote+=1
                        read_ids[read.query_name] = 1
                    i+=1
                    if i % 10000 == 0:
                        print(f"processed {i} reads in {time.time()-t} seconds and wrote {wrote} records. Currently at {chrom} {pos}")
                        t = time.time()
    outbam.close()
    bam.close()
    pysam.sort("-o", output_bam, tmp_bam)
    os.system(f"rm {tmp_bam}")
    pysam.index(output_bam)
    print("done")

def get_featuremap(output_bam, output_featuremap, tumor_vcf):
    output_featuremap_tmp = output_featuremap.replace('.vcf', '.tmp.vcf')
    
    cmd = f'bash workflow/scripts/WgsNoDuplexFeatureMap.sh {output_bam} {output_featuremap_tmp}'
    os.system(cmd)
    os.system(f'bedtools intersect -header -a {output_featuremap_tmp} -b {tumor_vcf} > {output_featuremap}')

    
def estimate_TF(output_featuremap, output_tf, sample_id, tumor_vcf):
    chromosomes = autosomal = ['chr'+str(i) for i in range(1,23)] #+ ['chrX', 'chrY']
    
    vcf = pd.read_csv(tumor_vcf, comment = '#', sep = '\t', names=['chrom', 'pos', 'id', 'ref', 'alt', 'a', 'b', 'c', 'd', 'e', 'f'])
    vcf = vcf[vcf['chrom'].isin(chromosomes)]

    df = pd.read_csv(output_featuremap, sep = '\t', names= ['chrom', 'pos', 'id', 'ref', 'alt', 'a', 'b', 'INFO'], comment='#')
    df['X_FILTEREDCOUNT'] = df['INFO'].str.split(';', expand=True)[4]
    df['X_FILTEREDCOUNT'] = df['X_FILTEREDCOUNT'].str.split('=', expand=True)[1].astype('int')
    
    df = df[['chrom', 'pos', 'ref', 'alt', 'X_FILTEREDCOUNT']]
    df2 = df.value_counts().reset_index()
    df2.columns = ['chrom', 'pos', 'ref', 'alt', 'X_FILTEREDCOUNT', 'alt_count']

    df3 = pd.merge(df2, vcf, on = ['chrom', 'pos', 'ref', 'alt'])
    
    # ctDNA
    ctDNA = sum(df3['alt_count']) / ( len(vcf.index) * np.mean(df3['X_FILTEREDCOUNT']) )
    
    # error
    error = len(df3[df3['alt_count']==1].index) / ( len(vcf.index) * np.mean(df3['X_FILTEREDCOUNT']) )
    
    # positions with a variant
    positions = len(df3[df3['alt_count']>=1].index)
    
    # vcf size
    size = len(vcf.index)
     
    with open(output_tf, 'w') as w:
        w.write('sample_id\ttumor_fraction\terror_rate\tpositions_with_alt\tvcf_size\n')
        w.write(f'{sample_id}\t{str(ctDNA)}\t{str(error)}\t{positions}\t{size}\n')


    
def main():

    input_bam, output_bam, output_featuremap, tumor_vcf, output_tf, sample_id = parse_arguments()
    intersect_tumor_bam(input_bam, output_bam, tumor_vcf)
    get_featuremap(output_bam, output_featuremap, tumor_vcf)
    
    estimate_TF(output_featuremap, output_tf, sample_id, tumor_vcf)
    
if __name__ == "__main__":
    main()
    
