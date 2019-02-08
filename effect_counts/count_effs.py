#!/usr/bin/env python

import sys
import re
import numpy as np

vcffile = sys.argv[1]

with open(vcffile, 'r') as vcf_handle:
    for line in vcf_handle:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                samples = {}
                sample_names = line.strip().split('\t')[9:]
                for i in sample_names:
                    samples['HIGH'] = [0] * len(sample_names)
                    samples['MODERATE'] = [0] * len(sample_names)
                    samples['LOW'] = [0] * len(sample_names)
                    samples['MODIFIER'] = [0] * len(sample_names)
#                print samples
            continue
        content = line.strip().split('\t')
        if '|HIGH|' in content[7]:
            for j in range(9, len(content)):
                if re.match('[01]/1', content[j][:3]):
                    samples['HIGH'][j - 9] += 1
            continue
        if '|MODERATE|' in content[7]:
            for j in range(9, len(content)):
                if re.match('[01]/1', content[j][:3]):
                    samples['MODERATE'][j - 9] += 1
            continue
        if '|LOW|' in content[7]:
            for j in range(9, len(content)):
                if re.match('[01]/1', content[j][:3]):
                    samples['LOW'][j - 9] += 1
            continue
        if '|MODIFIER|' in content[7]:
            for j in range(9, len(content)):
#                print len(samples['MODIFIER'])
#                print j
                if re.match('[01]/1', content[j][:3]):
                    samples['MODIFIER'][j - 9] += 1
            continue

print 'SAMPLE' + '\t' + 'TYPE' + '\t' + 'COUNT'
for i in range(len(sample_names)):
     for key in samples:
         print sample_names[i] + '\t' + key + '\t' + str(samples[key][i])
