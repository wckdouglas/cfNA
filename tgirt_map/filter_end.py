#!/usr/bin/env python

import pysam
import argparse
import six
from sequencing_tools.bam_tools import fragment_ends

def paired_bam(bam_handle):

    try:
        while True:
            #try:
                read1 = six.next(bam_handle)
                read2 = six.next(bam_handle)

                #assert(read1.is_read1 and read2.is_read2)
                assert(read1.query_name.split('/')[0] == read2.query_name.split('/')[0])
                yield read1, read2

    except StopIteration:
        pass


class filter_process:
    def __init__(self, inbam, tsRNA_type):
        self.inbam = inbam
        self.tsRNA_type = tsRNA_type
        self.small_RNA = 45
        self.buffer = 5
        self.ref_length = {ref:ref_len for ref_len, ref in zip(self.inbam.lengths,self.inbam.references)}

    def __five_prime_alignments(self):
        for aln1, aln2 in paired_bam(self.inbam):
            mapped = not (aln1.is_unmapped and aln2.is_unmapped)
            aln_start, aln_end = fragment_ends(aln1, aln2)
            five_prime = aln_start < self.buffer
            template_length = aln_end - aln_start
            small_RNA = template_length < self.small_RNA
            if mapped  and five_prime and small_RNA:
                yield aln1, aln2

    def __three_prime_alignments(self):
        for aln1, aln2 in paired_bam(self.inbam):
            rna = aln1.reference_name
            rna_end = self.ref_length[rna]
            aln_start, aln_end = fragment_ends(aln1, aln2)
            template_length = aln_end - aln_start

            three_prime = aln_end > (rna_end - self.buffer)
            mapped = not (aln1.is_unmapped and aln2.is_unmapped)
            small_RNA = template_length < self.small_RNA
            if mapped and three_prime and small_RNA:
                yield aln1, aln2


    def __fulllength_alignments(self):
        for aln1, aln2 in paired_bam(self.inbam):
            rna = aln1.reference_name
            rna_end = self.ref_length[rna]
            aln_start, aln_end = fragment_ends(aln1, aln2)
            template_length = aln_end - aln_start

            three_prime = aln_end > (rna_end - self.buffer)
            mapped = not (aln1.is_unmapped and aln2.is_unmapped)
            five_prime = aln_start < self.buffer
            if mapped and three_prime and five_prime:
                yield aln1, aln2

    def filter_alignments(self):
        if self.tsRNA_type == 'five':
            res = self.__five_prime_alignments()

        elif self.tsRNA_type == 'three':
            res = self.__three_prime_alignments()

        else:
            res = self.__fulllength_alignments()

        for aln1, aln2 in res:
            yield aln1, aln2


def getopt():
    parser = argparse.ArgumentParser(description = "Filter alignments with 5' or 3' tsRNA fragment")
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file (name-sorted; can be stdin, use -)')
    parser.add_argument('-o','--outbam', default = '-', help = 'output bam file (default: - )')
    parser.add_argument('--type', choices = ['five','three'],
                        help ='tsRNA type')
    args = parser.parse_args()
    return args

def main():
    args = getopt()
    with pysam.Samfile(args.inbam,'rb') as bam:
        with pysam.Samfile(args.outbam, 'wb', template = bam) as out:
            fp = filter_process(bam, args.type)
            for aln1, aln2 in fp.filter_alignments():
#                if aln1.reference_name.startswith('TR|MT-'):
                    out.write(aln1)
                    out.write(aln2)
                        

if __name__ == '__main__':
    main()
