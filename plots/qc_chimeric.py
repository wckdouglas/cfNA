#!/usr/bin/env python

import pysam
import numpy as np
import pandas as pd
import glob
from multiprocessing import Pool
import re
import os


forward_false_clip = re.compile('^[0][ACTG]')
reverse_false_clip = re.compile('[ACTG][0]$')


def count_softclipped(sample_folder):
    bam_file = sample_folder + '/Combined/primary.deduplicated.bam'
    samplename = os.path.basename(sample_folder)
    print('Running: %s' %bam_file)

    softclip = re.compile('[0-9]+S')
    softclip_count = 0
    bases = 0
    clipped_bases = 0
    aln_count = 0
    with pysam.Samfile(bam_file, 'r') as bam:
        for aln in bam:
            if not aln.is_duplicate:
                aln_count += 1
                bases += len(aln.query_sequence)
                if 'S' in aln.cigarstring:
                    softclip_count += 1
                    clips = softclip.findall(aln.cigarstring)
                    for clip in clips:
                        clipped_bases += int(clip.rstrip('S'))
    return (samplename, aln_count, softclip_count, bases, clipped_bases)


def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    out_path = project_path + '/softclip'
    out_table = out_path + '/clip_table.tsv'
    sample_folders = glob.glob(project_path + '/Q*001')
    sample_folders = filter(lambda x: not re.search('L[12]',x), sample_folders)
    p = Pool(24)
    rows = p.map(count_softclipped, sample_folders)
    p.close()
    p.join()

    df = pd.DataFrame(rows, columns = ['samplename','aln_count',
                                'softclip_count',
                                'bases','clipped_bases'])
    df.to_csv(out_table, sep='\t',index=False)
    print('Written %s' %out_table)
    




if __name__ == '__main__':
    main()
