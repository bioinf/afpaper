#!/usr/bin/env python3.6

import re, os, sys, getopt
from collections import Counter
from scipy.stats import binom
import codecs
import argparse
from collections import defaultdict

def main(args):
    required_info_fields = ("AF,AN,ANN,"
                            "AF_afr,AN_afr,AF_sas,AN_sas,AF_amr,AN_amr,AF_eas,AN_eas,AF_nfe,AN_nfe,AF_fin,AN_fin,"
                            "GNOMAD_AF,CLNSIG,CLNDN")
    required_info_fields = required_info_fields.split(',')
    binom_subset = [i for i in required_info_fields if all([re.search(f'AF', i)!=None, 'AF'!=i])]

    df_dict = defaultdict(list)
    with codecs.open(args.input, encoding='utf-8', errors='ignore') as f:
        for line in f:
            # getting header
            if '#CHROM' == line[:6]:
                sample_names = line.rstrip().split()[9:]
                break
        # skip header and iterate body
        for line in f:
            if '#' not in line:
                line = line.rstrip().split()

                info = line[7]

                df_dict['chrom'].append(line[0])
                df_dict['pos'].append(int(line[1]))
                df_dict['rsID'].append(line[2])
                df_dict['Ref'].append(line[3])
                df_dict['Alt'].append(line[4])
                for i in required_info_fields:    
                    try:
                        df_dict[i].append(re.search(f"{i}=([^;]+)", info).group(1))
                    except:
                        df_dict[i].append('.')

                gt = line[9:]
                # cat genotype info
                for index, item in enumerate(gt):
                    gt[index] = gt[index][:3]


                # counting genotypes
                gt_dict = Counter(gt)
                try:
                    del gt_dict['./.']
                except KeyError:
                    pass

                # getting AN
                an = sum(gt_dict.values()) * 2

                # getting indices
                indicies_1_1 = [i for i, s in enumerate(gt) if '1/1' == s[:3]]
                indicies_0_1 = [i for i, s in enumerate(gt) if '0/1' == s[:3]]

                # getting names of samples with particular gt
                names_samples1_1, names_samples0_1 = [], []

                for i in indicies_0_1:
                    names_samples0_1.append(sample_names[i])

                for j in indicies_1_1:
                    names_samples1_1.append(sample_names[j])

                names_samples1_1, names_samples0_1 = "|".join(names_samples1_1), "|".join(names_samples0_1)
                df_dict['HOM'].append(gt_dict["1/1"])
                df_dict['NAMES_HOM_SAMPLES'].append(names_samples1_1)
                df_dict['HET'].append(gt_dict['0/1'])
                df_dict['NAMES_HET_SAMPLES'].append(names_samples0_1)
                # get binom/poisson p-value

                for af in binom_subset:
                    p_name = af.strip('AF').strip('_')
                    try:
                        num_alleles = len(indicies_0_1) + 2 * len(indicies_1_1)

                        binom_pval = binom.pmf(num_alleles, int(df_dict['AN'][-1]), float(df_dict[af][-1]))
                        df_dict[f'{p_name.upper()}_P_biom'].append(binom_pval)
                    except:
                        df_dict[f'{p_name.upper()}_P_biom'].append(1)



    df = pd.DataFrame(df_dict)
    df['gene'] = df['ANN'].str.split('|').str[3]

    df['GNOMAD_AF'] = df['GNOMAD_AF'].convert_objects(convert_numeric=True)
    df['AF'] = df['AF'].convert_objects(convert_numeric=True)

    known_pathogenic = df[(df['CLNSIG'].str.contains('Pathogenic'))&(df['GNOMAD_AF']<0.005)]
    expected_pathogenic = df[(df['ANN'].str.contains('stop_gained|frameshift|splice_donor_variant|splice_acceptor_variant'))&
       (df['AF']<0.01)]

    known_pathogenic.to_csv('known_pathogenic.tsv', sep='\t', index=None)
    expected_pathogenic.to_csv('expected_pathogenic.tsv', sep='\t', index=None)
    

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input', type=str,
                        help='Path to input vcf file.')

    args = parser.parse_args()
    main(args)

