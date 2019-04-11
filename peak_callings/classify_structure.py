#!/usr/bin/env python

import RNA
from skbio import io
import re
import sys

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <fasta file>' %sys.argv[0])
fa = sys.argv[1]
for r in io.read(fa, 'fasta'):
    seq = str(r)[20:-20]
    f, e = RNA.fold(seq.strip('N'))
    folded = RNA.b2C(f)
    is_cloverleave = re.findall('[A-Z]', folded)
    is_tRNA = is_cloverleave and 'HHH' in ''.join(is_cloverleave) 
    closed_end = folded.startswith('(') and folded.endswith(')')
    cloverleave = 'cloverleaf' if is_tRNA and closed_end else 'hairpin'
    print(r.metadata['id'], cloverleave, folded, seq.strip('N'))


        
