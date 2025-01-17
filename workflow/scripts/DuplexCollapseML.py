#!/usr/bin/env python
# coding: utf-8


import pysam
import pandas as pd
import argparse
import re
from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn import metrics
import numpy as np
import pickle
pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_columns', None)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Add a couple features you know how it is')
    parser.add_argument('--AnnotatedFeatureMap', help = 'AnnotatedFeatureMap file with the necessary tags')
    parser.add_argument('--output')
    parser.add_argument('--model')
    args = parser.parse_args()
    featuremap_file = args.AnnotatedFeatureMap
    output_file = args.output
    model = args.model
    
    return(featuremap_file, output_file, model)


def DuplexCollapse(featuremap_file, output_file, model):
    clf = pickle.load(open(model, 'rb'))
    df = pd.read_csv(featuremap_file, sep = '\t', na_filter=False)
    print(len(df.index))
    #print(df)
    df = df[df['MI'] != -1]
    #print(df)
    print(len(df.index))
    # Easy case where FM_passed == cU
    print(df[['FM_passed', 'cU']])
    df[['FM_passed', 'cU', 'X_RN']].to_csv('test')
    df2 = df[df['FM_passed'] == df['cU']]
    print(len(df2.index))
    df2['X_FLAGS'].replace({1024: 0, 1040: 1, 16: 1}, inplace=True)

    # duplexable reads
    duplex_MI = df2.loc[df2['X_FLAGS'] == 0, 'MI'][df2.loc[df2['X_FLAGS'] == 0, 'MI'].isin(df2.loc[df2['X_FLAGS'] == 1, 'MI'])].unique()
    df3 = df2[df2['MI'].isin(duplex_MI)]

    df3 = df3[["alt", "X_FLAGS", "duplex_bp", "MI_pos"]]

    for base in ['A', 'T', 'C', 'G', 'N']:
        df3[f'{base}_top'] = (df3['X_FLAGS'] == 0) & (df3['alt'] == base)
        df3[f'{base}_bot'] = (df3['X_FLAGS'] == 1) & (df3['alt'] == base)

    df3.drop(columns='alt', inplace=True)

    df4 = df3.groupby(['MI_pos', 'duplex_bp']).sum().reset_index()

    features = df4[['A_top', 'T_top', 'C_top', 'G_top', 'N_top', 'A_bot', 'T_bot', 'C_bot', 'G_bot', 'N_bot']]
    
    y_pred = (clf.predict(features))
    df4['duplex_bp2'] = y_pred
    df4.drop(columns = 'X_FLAGS', inplace=True)
    df_final = pd.merge(df, df4, on = ['MI_pos', 'duplex_bp'])
    
    df_final.to_csv(output_file, sep = '\t', index = False)

def main():
    featuremap_file, output_file, model = parse_arguments()    
    DuplexCollapse(featuremap_file, output_file, model)



if __name__ == '__main__':
    main()

