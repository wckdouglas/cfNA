#!/usr/bin/env python

import pandas as pd
import os
import sys
import numpy as np
import fileinput
from operator import itemgetter

line_template = '{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\tRBP\t{name}'
for line in fileinput.input():
    fields = line.split('\t')
    chrom, start, end, name, strand, pval = itemgetter(0,1,2,3,5,7)(fields)
    try:
        rbp, cell, rep = itemgetter(0,-2,-1)(name.split('_'))
    except ValueError:
        sys.exit(name)
    if float(pval) > 2:
        print(line_template.format(chrom = chrom,
                             start = start,
                             end = end,
                             name = rbp,
                             score = pval,
                             strand = strand))

