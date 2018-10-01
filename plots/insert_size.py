#!/bin/env python

import glob
import os
import pandas as pd
from collections import defaultdict, Counter
from functools import reduce
from multiprocessing import Pool
import re
from sequencing_tools.io_tools import xopen
from functools import partial

def add_dict(a, b):
    keys = set(a.keys()).union(b.keys())
    for key in keys:
        a[key] += b[key]
    return a


def bed_fragments(bed_file):
    size_dict = Counter()
    
    with xopen(bed_file, 'r') as bed:
        for line in bed:
            fields = line.split('\t')
            isize = int(fields[2]) - int(fields[1])
            size_dict[isize] += 1
    return size_dict

def sample_fragments(sample_folder):
    tRNA_bed_file = sample_folder + '/smallRNA/aligned.dedup.bed.gz'
    rRNA_bed_file = sample_folder + '/rRNA_mt/aligned.dedup.bed.gz'
    samplename = os.path.basename(sample_folder)


    print('Running %s' %samplename)
    bed_folder = os.path.dirname(sample_folder) + '/bed_files'

    all_bed_file = bed_folder + '/' + samplename + '.bed.gz'

    size_dict = bed_fragments(tRNA_bed_file)
    size_dict += bed_fragments(rRNA_bed_file)
    size_dict += bed_fragments(all_bed_file)
    return size_dict

def get_isize(samples, insert_size_path, args): 
    regex, label = args
    samplenames = filter(lambda x: re.search(regex, os.path.basename(x)), samples)
    size_dicts = map(sample_fragments, samplenames)
    size_dict = reduce(add_dict, size_dicts)

    out_file = insert_size_path + '/' + label + '.tsv'
    pd.DataFrame({'isize': list(size_dict.keys()),
                  'size_count': list(size_dict.values())})\
        .to_csv(out_file, index=False, sep='\t')
    print('Written: ', out_file)



def main():
    project_path='/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    insert_size_path = project_path + '/fragment_sizes'
    if not os.path.isdir(insert_size_path):
        os.mkdir(insert_size_path)
    samples = glob.glob(project_path + '/*001')

    isize_func = partial(get_isize, samples, insert_size_path)

    args = zip(['Q[Cc][Ff][0-9]+_', 'Frag[1123]', '_L[12]_', 'N[aA][1-9]+', 
                             '[aA]ll[12]', '[DE][ED]|Exo','_HS[0-9]+_'],
                            ['unfragmented','fragmented','polyA','alkaline', 
                             'all', 'Exo','high_salt'])
    p = Pool(24)
    p.map(isize_func, args)
    p.close()
    p.join()
    return 0

if __name__ == '__main__':
    main()

        




    





