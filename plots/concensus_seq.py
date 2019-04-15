#!/usr.bin/env python

import subprocess
import numpy as np
import pyBigWig as pbw

class concensus():
    def __init__(self, bam_file, coverage_files = [], debug=False):
        '''
        testing:
        workdir = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
        fwd = workdir + '/bed_files/merged_bed/coverage/unfragmented.fwd.bigWig'
        rvs = workdir + '/bed_files/merged_bed/coverage/unfragmented.rvs.bigWig'
        bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam'
        con = concensus(bam, coverage_files = [fwd, rvs])
        '''
        self.bam = bam_file
        self.debug = debug
        self.fwd_coverage = pbw.open(coverage_files[0])
        self.rvs_coverage = pbw.open(coverage_files[1])

    def __extract_concensus__(self, coor, strand):
        '''
        extract all Read2 sequences from peak summit and multiple-aligned,
        find concensus
        '''
        strand_flag = '0x20' if strand == '-' else '0x10'
        command = "samtools view -b {bam} '{coor}'"\
            '| samtools sort -n | samtools fastq -F0x40 -f {flag} - '\
            '| seqtk seq -a | muscle  | seqtk seq '\
            "| python muscle_to_concensus.py | rev | tr 'ACTGU' 'TGACA' "\
                .format(bam = self.bam, coor = coor, flag = strand_flag)
        if self.debug:
            print(command)
        ps = subprocess.Popen(command, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        res =  ps.communicate()[0]
        seq =  res.decode().split('\n')[-2]
        seq = '' if 'ValueError' in seq else seq
        return seq


    
    def __narrow_coor__(self, chrom, start, end, strand):
        '''
        find peak summit
        '''
        bw = self.fwd_coverage if strand == '+' else self.rvs_coverage
        coverage = bw.values(chrom, start, end, numpy =True)
        start += coverage.argmax()
        return chrom + ':' + str(start) + '-' + str(start + 1)
        
    def find_concensus(self, chrom, start, end ,strand):
        coor = self.__narrow_coor__(chrom, start, end, strand)
        return self.__extract_concensus__(coor, strand)


