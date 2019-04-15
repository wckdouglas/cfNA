#!/usr/bin/env python

import pandas as pd
import numpy as np
from functools import partial
import re
import io


def parse_line(field_starts, field_ends, line, header = False):
    if header:
        return [line.strip('\n#')[s:e].strip() for s, e in zip(field_starts, field_ends)]
    else:
        if len(line.split(' ')[0]) < field_ends[0] - field_starts[0]:
            return [line.strip('\n#')[s:e+1].strip() for s, e in zip(field_starts, field_ends)]


def define_table(infile):
    header = next(infile)
    field_width = next(infile)
    field_widths = field_width.strip('#\n').split(' ')
    field_widths = np.cumsum(np.array(list(map(len, field_widths))) + 1)
    field_starts = np.roll(np.append(field_widths,[0]),1)
    header = parse_line(field_starts, field_widths, header, header=True)
    return header, field_starts, field_widths




def read_tbl(tbl_file):
    with open(tbl_file) as infile:
        records = []
        header, field_starts, field_ends = define_table(infile)
        line_parser = partial(parse_line, field_starts, field_ends)
        for line_count, line in enumerate(infile):
            if not line.startswith('#'):
                fields = line_parser(line)
                if fields:
                    record = {h:f for f, h in zip(fields, header)}
                    records.append(record)
    return pd.DataFrame(records) \
            .assign(score = lambda d: d.score.str.extract('([0-9]+\.[0-9]+)$', expand=False)) \
            .assign(score=lambda d: d.score.fillna('0').astype(float))


def read_tbl_old(tbl_file, truncate=None):
    with open(tbl_file) as infile:
        header = 'target name|accession|query name|accession_null|mdl|mdl from|mdl to'\
                '|seq from|seq to|strand|trunc|pass|'\
                'gc|bias|score|E-value|inc|description of target'
        headers = header.split('|')
        rows = []
        for i, line in enumerate(infile):
            if not line.startswith('#'):
                fields = line.strip().split(' ')
                fields = list(filter(lambda x: x!='',fields))
                rows.append(fields)
            if truncate and i ==truncate:
                break
    return pd.DataFrame(rows, columns = headers)
        
                

def read_tbl(tbl_file):
    lines = []
    header = ['target name',
            'accession',
            'query name',
            'accession strand',
            'mdl',
            'mdl from',   
            'mdl to',
            'seq from',
            'seq to',
            'strand',
            'trunc',
            'pass', 'gc','bias', 
            'score','E-value','inc', 'description of target']
    with open(tbl_file) as infile:
        for line in infile:
            if not line.startswith('#'):
                line = line.strip('#').strip()
                lines.append(re.sub('\s+','\t', line))
    return pd.read_table(io.StringIO('\n'.join(lines)), names = header)


                


