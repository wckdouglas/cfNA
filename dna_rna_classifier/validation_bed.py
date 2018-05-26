#!/usr/bin/env python

import sys
from collections import defaultdict
import random

if len(sys.argv) != 4:
    sys.exit('[usage] python %s <inbed> <out_prefix> <test_number>' %(sys.argv[0]))


out_prefix = sys.argv[2]
test_bed = out_prefix + '/test.bed'
train_bed = out_prefix + '/train.bed'
label_counter = defaultdict[int]
test_number = int(sys.argv[3])
half_test = test_number//2

test_count = 0
train_count = 0
with open(sys.argv[1],'r') as inbed,\
        open(test_bed, 'w') as test,\
        open(train_bed, 'w') as train:
    for line_num, bed_line in enumerate(inbed):
        bed_line = bed_line.strip()
        fields = bed_line.split('\t')

        if 'SRR' not in fields[3] \
               and label_counter[fields[-1]] < half_test \
            print(bed_line, file = test)
            test_count += 1

        else:
            print(bed_line, file = train)
            train_count += 1
message  = 'Parsed {all} lines\n'\
        'Test sample: {test}\n'\
        'Train sample: {train}'.format(all = line_num, test = test_count, train = train_count) 
