#!/usr/bin/env python

import os

project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup'
coors = ['chr22:24,349,641-24,349,711', 'chr10:71,355,032-71,355,103',
        'chr1:181,392,082-181,392,121', 'chr10:101817588-101817660','chr4:156,379,949-156,380,020']
names = ['chr22_unkonwn', 'chr10_unknown',
        'CACNA1E', 'CPN1','chr4_unkown']
for coor, name in zip(coors, names):
    command = "samtools view -b {path}/unfragmented.chrM_filter.dedup.bam '{coor}'"\
            '| samtools sort -n | samtools fastq  - '\
            '| pe_fq_merge.py -a -e 0.2 '\
            '| seqtk seq -a | muscle  | seqtk seq '\
            '| python muscle_to_concensus.py '\
            "| RNAshapes -s -e 10 -l -m \'[[][][]]\' "\
            '> {path}/fold/{name}.shapes'\
            .format(path = project_path, coor = coor, name = name)
    print(command)
#    os.system(command)

