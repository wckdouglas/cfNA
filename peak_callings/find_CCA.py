#!/usr/bin/env python

import os

project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup'
coors = ['chr22:24,349,641-24,349,711', 'chr10:71,355,032-71,355,103',
        'chr1:181,392,082-181,392,121', 'chr10:101817588-101817660',
        'chr9:5095156-5095227',
        'chr4:156,379,949-156,380,020',
        'chr2:202,077,805-202,077,873']
names = ['chr22_unkonwn', 'chr10_unknown',
        'CACNA1E', 'CPN1',
        'JAK2',
        'chr4_unkown',
        'CASP10']


files = []
for coor, name in zip(coors, names):
    outfile =  project_path + '/fold/' + name + '.shapes'
    fa = outfile.replace('.shapes','.fa')
    command = "samtools view -b {path}/unfragmented.chrM_filter.dedup.bam '{coor}'"\
            '| samtools sort -n | samtools fastq -F0x40  - '\
            '| seqtk seq -a | muscle  | seqtk seq '\
            '| python muscle_to_concensus.py '\
            "| rev | tr 'ACTGU' 'TGACA' " \
            "| RNAshapes -s -e 10 -l -m \'[[][][]]\' "\
            '> {path}/fold/{name}.shapes'\
            .format(path = project_path, coor = coor, name = name)
     
    print(command)
    os.system(command)
    with open(outfile) as inf,\
            open(fa,'w') as fasta:
        seq = next(inf).strip()
        print('>{name}\n{seq}'.format(seq = seq.replace('T','U'), name = name), file = fasta)
    files.append(fa)


combined_fa = project_path + '/fold/unknown_tRNA.fa'
command = 'cat {files} > {out} '.format(files = ' '.join(files),
                                         out = combined_fa)
os.system(command)







