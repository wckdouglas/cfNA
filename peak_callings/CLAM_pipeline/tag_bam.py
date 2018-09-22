#!/usr/bin/env python

from itertools import groupby
import sys
import pysam
import argparse
from tqdm import tqdm

def getopt():
    parser = argparse.ArgumentParser(description = 'adding tag to bam')
    parser.add_argument('-i', '--inbam', required=True, help = 'input bam file')
    parser.add_argument('-o','--outbam', default = '-', help = 'output bam file')
    parser.add_argument('-t','--tag',default='NH', help = 'tag name')
    parser.add_argument('--tag_value', default='', help = 'tag value to add (required if not tag != NH)')
    args = parser.parse_args()
    return args


def process_NH(inbam, outbam):
    for read_name, alns in groupby(inbam, lambda a: a.query_name):
        alns = list(alns)

        NH = int(len(alns)/2)
        for aln in alns:
            if not aln.is_unmapped:
                aln.tags += [('NH',NH)]
        outbam.write(aln)



def process_TAG(inbam, outbam, tag, tag_value):
    for aln in tqdm(inbam):
        aln.set_tag(tag, tag_value, replace=True)
        outbam.write(aln)


def main():
    args = getopt()
    assert args.tag_value or args.tag == "NH", 'Eiher provide a tag_value or use tag = NH (number of hits)'

    with pysam.Samfile(args.inbam, 'rb') as inbam:
        with pysam.Samfile(args.outbam,'wb', template=inbam) as outbam:
            if args.tag == "NH":
                process_NH(inbam, outbam)

            else:
                process_TAG(inbam, outbam, args.tag, args.tag_value)

if __name__ == '__main__':
    main()


