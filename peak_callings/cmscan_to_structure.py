#!/usr/bin/env python

import fileinput
import os
import re

seq_re = re.compile('[ACTGUactgu*\[\]0-9\>\<-]+')
def get_seq(line):
    res = seq_re.findall(line)
    longest = max(map(len, res))
    return list(filter(lambda x: len(x)==longest, res))[0]



filename = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated/unfragmented.Long_RNA.cmscan'
rfam=''
with open(filename) as f:
    while True:
        try:
            line = next(f)
            if line.startswith('Query'):
                peakname = line.strip().split('       ')[1].split(' ')[0]

            elif line.startswith('>>'):
                if rfam != '' and '[' not in seq:
                    seq = re.sub('[0-9]+', '',seq)
                    fold = re.sub('[~]+$','', fold)
                    assert len(seq) == len(fold), seq +'\t' + fold + '\t%i,%i' %(len(seq),len(fold))
                    print(peakname, rfam, seq, fold)
                
                rfam = line.strip().split(' ')[1]
                fold = ''
                ref = ''
                seq = ''

            elif line.strip().endswith('CS'):
                fold += line.strip().replace(' CS','')
                ref += get_seq(next(f).strip())
                line = next(f).strip()
                seq += get_seq(next(f).strip())
        except StopIteration:
            print(peakname, rfam, seq, fold)

        





