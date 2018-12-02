#!/usr/bin/env python

import pysam
import numpy as np
import pandas as pd
import glob
from multiprocessing import Pool
import re
import os
from collections import defaultdict
from sequencing_tools.fastq_tools import reverse_complement


def extract_non_template(aln):
    cut_pos, cut_base, cut_qual = [], '',[]
    lpos = -100

    for (seq_pos, ref_pos, ref_base) in aln.get_aligned_pairs(with_seq=True):
        if seq_pos and (not ref_base or ref_base != aln.query_sequence[seq_pos]):
            base, qual = aln.query_sequence[seq_pos], aln.query_qualities[seq_pos]

            if  lpos>=0 and seq_pos != lpos + 1:
                yield (cut_pos, cut_base, cut_qual)
                cut_pos = []
                cut_base = ''
                cut_qual = []

            cut_pos += [seq_pos]
            cut_base += base
            cut_qual += [qual]
        

            lpos = seq_pos

    if cut_pos:
        yield (cut_pos, cut_base, cut_qual)

def count_softclipped(sample_folder):
    '''
    Sense strand:
     R1: S------> 
         S<------- :R2

    Antisense strand:
     R2: ------>S 
         <------S :R1


    '''
    # sample_folder = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/Qcf10_R1_001'
    bam_file = sample_folder + '/Combined/primary.duplicated.bam'
    samplename = os.path.basename(sample_folder)
    print('Running: %s' %bam_file)

    non_template_additions = defaultdict(int)
    with pysam.Samfile(bam_file, 'r') as bam:
        for aln in bam:
            if not aln.is_duplicate and aln.is_read1:
                non_template_added = 0

                for seq_poses, bases, quals in extract_non_template(aln):
                    if aln.is_reverse and max(seq_poses) == len(aln.query_alignment_sequence): 

                        #antisense
                        #if min(quals) > 20:
                            non_template_additions[reverse_complement(bases)] += 1
                            non_template_added = 1

                    elif not aln.is_reverse and min(seq_poses) == 0:
                        #sense
                        #if min(quals) > 20:
                            non_template_additions[bases] += 1
                            non_template_added = 1
                    
                if not non_template_added: 
                    non_template_additions['.'] += 1
    

    return pd.DataFrame({'non_template': list(non_template_additions.keys()),
                'non_template_count': list(non_template_additions.values()),
                'samplename': samplename})



def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    out_path = project_path + '/non_template_added'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    out_table = out_path + '/non_template_table.feather'
    sample_folders = glob.glob(project_path + '/Q*001')
    sample_folders = filter(lambda x: not re.search('L[12]',x), sample_folders)
    p = Pool(24)
    rows = p.map(count_softclipped, sample_folders)
    p.close()
    p.join()

    pd.concat(rows).to_feather(out_table)
    print('Written %s' %out_table)
    




if __name__ == '__main__':
    main()
