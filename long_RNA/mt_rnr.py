#!/usr/bin/env python

from pybedtools import BedTool
import mappy
import os
from itertools import groupby
from sequencing_tools.fastq_tools import reverse_complement


genes_bed = BedTool(os.environ['REF'] + '/hg19/new_genes/genes.bed') \
        .filter(lambda x: x[3] in ['MTRNR2L12','MTRNR2L8'])


bam_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bam'
bam = BedTool(bam_path + '/unfragmented.name_sort.bam') \
    .pair_to_bed(b = genes_bed, type="either") \


mt_count = 0
chrM = mappy.Aligner(os.environ['REF'] + '/hg19/genome/chrM.minimap2_idx',
                     preset = 'sr')

out_sam = 'mtrnr_mt.sam'
sam = open(out_sam,'w')
for aln_count, (seq_name, alns) in enumerate(groupby(bam, key=lambda x: x[0])):
    aln1, aln2 = alns

    if int(aln1[1]) & 0x40:
        aln1, aln2 = aln1, aln2  
    else: 
        aln1, aln2 = aln2, aln1

    seq1 = reverse_complement(aln1[9]) if int(aln1[1])&0x16 else aln1[9]
    seq2 = reverse_complement(aln2[9]) if int(aln2[1])&0x16 else aln2[9]


    aln = chrM.map(seq1, seq2)

    try:
        a = next(aln)
        mt_count += 1
        aln1[9],aln2[9] = seq1, seq2
        sam.write(str(aln1) +str(aln2))
    except StopIteration:
        pass

print('MT: %i/%i reads' %(mt_count,aln_count))
sam.close()






    

