#!/usr/bin/env python

import sys
from collections import defaultdict
import random

if len(sys.argv) != 4:
    sys.exit('[usage] python %s <inbed> <out_prefix> <test_number>' %(sys.argv[0]))


out_prefix = sys.argv[2]
test_bed = out_prefix + '/test.bed'
train_bed = out_prefix + '/train.bed'
train_RNA_bed = out_prefix + '/train_RNA.bed'
train_DNA_bed = out_prefix + '/train_DNA.bed'
label_counter = defaultdict(int)
test_number = int(sys.argv[3])
half_test = test_number//2

test_count = 0
train_count = 0
DNA = 0
RNA = 0
chroms = list(range(1,23))
chroms.append(['X','Y'])
chroms = ['chr{}'.format(c) for c in chroms]
with open(sys.argv[1],'r') as inbed,\
        open(test_bed, 'w') as test,\
        open(train_bed, 'w') as train,\
        open(train_RNA_bed, 'w') as train_RNA,\
        open(train_DNA_bed, 'w') as train_DNA:
    for line_num, bed_line in enumerate(inbed):
        bed_line = bed_line.strip()
        fields = bed_line.split('\t')
        chrom = fields[0]
        if chrom in chroms:
            label = fields[-1]

            if not fields[3].startswith('SRR') \
                   and label_counter[label] <= half_test:
                print(bed_line, file = test)
                label_counter[label] += 1
                test_count += 1

            else:
                print(bed_line, file = train)
                train_count += 1
                if bed_line.endswith('DNA'):
                    print(bed_line, file = train_DNA)
                    DNA += 1
                else:
                    print(bed_line, file = train_RNA)
                    RNA += 1
assert(DNA + RNA == train_count)
message  = 'Parsed {all} lines\n'\
        'Test sample: {test}\n'\
        'Train sample: {train}'.format(all = line_num, 
                                       test = test_count, 
                                       train = train_count) 
print(message)
