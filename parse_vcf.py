#!/usr/bin/env python3.6

import re, os, sys, getopt
from collections import Counter
from scipy.stats import binom
from math import ceil


def try_to_get_required_field(pattern, info):

    if pattern == 'GNOMAD_AF=([e\d\.\-]+)' or pattern == 'AF=([e\d\.\-]+)':
        try:
            field = re.search(pattern, info).group(1)
        except AttributeError:
            field = 0
    else:
        try:
            field = re.search(pattern, info).group(1)
        except AttributeError:
            field = "."

    return field

def main():

    # full_path = os.path.dirname(os.path.realpath(__file__))
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])

    # Process help
    for o, a in opts:
        if o in ("-h", "--help"):
            return 'NO HELP HERE'

    if len(args) != 2:
        print("ARGS: VCF_FILE OUTPUT")
        return

    vcf, out = args

    ### just for debug ###
    #  vcf, out = 'sample.tsv', '.'
    try:
        os.mkdir(f'{out}')
    except:
        pass


    with open(vcf, 'r+') as f, open(f'{out}/known_pathogenic_GNOMAD.tsv', 'w+') as kp, \
                               open(f'{out}/expected_pathogenic_GNOMAD.tsv', 'w+') as ep:

        # getting sample names
        for line in f:
            line = line.rstrip()
            if '#CHROM' == line[:6]:
                sample_names = line.split('\t')[9:]
                break
        # print(sample_names)

        for line in f:
            if '#' != line[:1]:

                line = line.rstrip()
                splitter_vcf = line.split('\t')

                # getting columns 10+ from .vcf
                gt = splitter_vcf[9:]
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

                # getting info field from .vcf
                info = splitter_vcf[7]

                # if it exist, we are taking required field
                ann = try_to_get_required_field('ANN=([^;]+)', info)
                clnsig = try_to_get_required_field('CLNSIG=([ \w\._]+)', info)
                clndn = try_to_get_required_field('CLNDN=([^;]+)', info)
                gnomad_score = try_to_get_required_field('GNOMAD_AF=([e\d\.\-]+)', info)
                af_score = try_to_get_required_field(';AF=([e\d\.\-]+)', info)

                try:
                    num_alleles = len(indicies_0_1) + 2 * len(indicies_1_1)
                    print(num_alleles)
                    binom_pval = binom.pmf(num_alleles, an, float(gnomad_score))
                except:
                    binom_pval = 1

                # print(binom_pval, an, af_score, gnomad_score)
                splitter_info = info.split(';')

                # getting known_pathogenic file
                var = re.split("\|", ann)

                if clnsig == 'Pathogenic':
                    try:
                        if float(gnomad_score) < 0.005:
                            kp.write(f'{splitter_vcf[0]}\t{splitter_vcf[1]}\t'
                                     f'{splitter_vcf[2]}\t{splitter_vcf[3]}\t'
                                     f'{splitter_vcf[4]}\t{gnomad_score}\t'
                                     f'{splitter_info[0][3:]}\t{binom_pval}\t{var[3]}\t'
                                     f'{gt_dict["0/1"]}\t{names_samples0_1}\t'
                                     f'{gt_dict["1/1"]}\t{names_samples1_1}\t'
                                     f'{an}\t{clndn}\t{clnsig}\t{var[1]}\n')
                    except ValueError:
                        kp.write(f'{splitter_vcf[0]}\t{splitter_vcf[1]}\t'
                                 f'{splitter_vcf[2]}\t{splitter_vcf[3]}\t'
                                 f'{splitter_vcf[4]}\t{gnomad_score}\t'
                                 f'{splitter_info[0][3:]}\t{binom_pval}\t{var[3]}\t'
                                 f'{gt_dict["0/1"]}\t{names_samples0_1}\t'
                                 f'{gt_dict["1/1"]}\t{names_samples1_1}\t'
                                 f'{an}\t{clndn}\t{clnsig}\t{var[1]}\n')
                # getting expected_pathogenic file
                if any(effect in ann for effect in ("stop_gained", "frameshift",
                                                    "splice_donor_variant",
                                                    "splice_acceptor_variant")):
                    if float(af_score) < 0.01:
                        ep.write(f'{splitter_vcf[0]}\t{splitter_vcf[1]}\t'
                                 f'{splitter_vcf[2]}\t{splitter_vcf[3]}\t'
                                 f'{splitter_vcf[4]}\t{gnomad_score}\t'
                                 f'{splitter_info[0][3:]}\t{binom_pval}\t{var[3]}\t'
                                 f'{gt_dict["0/1"]}\t{names_samples0_1}\t'
                                 f'{gt_dict["1/1"]}\t{names_samples1_1}\t'
                                 f'{an}\t{clndn}\t{clnsig}\t{var[1]}\n')

    col_names = ["Chr", "Pos", "rsID", "Ref", "Alt", "GNOMAD_AF", "AF",
                 "BINOM_p_value", "GENE_NAME", "HET", "NAMES_HET_SAMPLES",
                 "HOM", "NAMES_HOM_SAMPLES", "AN", "CLNDN", "CLNSIG", "ANN"]
    col_names = '\t'.join(map(str, col_names)) + "\n"

    with open(f'{out}/known_pathogenic_GNOMAD.tsv', 'r+') as kpr, \
            open(f'{out}/expected_pathogenic_GNOMAD.tsv', 'r+') as epr:
        k, e = kpr.read(), epr.read()
        with open(f'{out}/known_pathogenic_GNOMAD.tsv', "w+") as kpw, \
                open(f'{out}/expected_pathogenic_GNOMAD.tsv', "w+") as epw:
            kpw.write(f"{col_names}{k}")
            epw.write(f"{col_names}{e}")

if __name__ == '__main__':
    main()
