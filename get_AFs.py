import os
import re
import pandas as pd
import numpy as np
import gzip
import allel
from collections import defaultdict
import codecs
import argparse


def main(args):
    
    s = 'AF_afr,AN_afr,AF_sas,AN_sas,AF_amr,AN_amr,AF_eas,AN_eas,AF_nfe,AN_nfe,AF_fin,AN_fin,GNOMAD_AF'
    info_fields = s.split(',')
    AF_list = [i for i in info_fields if 'AF_' in i]
    
    df_dict = defaultdict(list)
    with codecs.open(args.vcf_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if '#' not in line:
                line = line.rstrip().split('\t')
                if all(['AF_' in line[7], 'GNOMAD_AF' in line[7]]):
                    info = line[7].split(';')
                    df_dict['chr'].append(line[0])
                    df_dict['pos'].append(line[1])
                    df_dict['rs'].append(line[2])
                    df_dict['AF'].append(info[1].split('=')[1])
                    d = {j.split('=')[0]: j.split('=')[1] for j in [i for i in info if any(['AF_' in i, 'GNOMAD_AF' in i])]}
                    d = {**d, **{i:'0.0' for i in set(AF_list).difference(set(d.keys()))}}
                    for k, v in d.items():
                        df_dict[k].append(v) 
                        
    df = pd.DataFrame(df_dict)
    df = df[~df['AF'].str.contains('\|')]
    
    for j in [i for i in df.columns if 'AF' in i] + ['pos']:
        df[j] = pd.to_numeric(df[j], errors='coerce')
        
    df.to_csv(args.output, sep='\t', index=None)
#     return df_dict
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-vcf_path', type=str,
                        help='Path to input vcf file (can be in .gz format).')
    parser.add_argument('-output', type=str, help='Path to output table.')

    args = parser.parse_args()
    main(args)
#     class Args:
        
#         def __init__(self):
#             self.vcf_path = '../data/gnomad_annotated.vcf'
#             self.output = ''
#     args = Args()
#    d = main(args)

