#!/usr/bin/env python

import glob
import os
import pandas as pd
from tblout_parser import read_tbl

project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated'
tbl_files = glob.glob(project_path + '/*tblout')
tbl = pd.concat(map(lambda x: read_tbl(x).assign(peak_type = lambda d: os.path.basename(x).split('.')[1]), tbl_files))
tbl.pipe(lambda d: d[d['E-value']<0.01]).to_csv('infernal.csv', index=False)


