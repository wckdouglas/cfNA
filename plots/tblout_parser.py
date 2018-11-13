#!/usr/bin/env python

import pandas as pd
import numpy as np
from functools import partial


def parse_line(field_starts, field_ends, line, header = False):
    if header:
        return [line.strip('\n#')[s:e].strip() for s, e in zip(field_starts, field_ends)]
    else:
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
                record = {h:f for f, h in zip(fields, header)}
                records.append(record)
    return pd.DataFrame(records) 



                

