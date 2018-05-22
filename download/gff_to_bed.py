#!/usr/bin/env python

import sys
from operator import itemgetter
from urllib.request import Request, urlopen
import gzip


def parse_extra_fields(extra_field):
    info_fields = extra_field.split(';')
    info_dict = {}
    for info in info_fields[:-1]:
        row_fields = info.strip().split(' ')
        info_dict[row_fields[0]] = row_fields[1].strip('')

    return info_dict


gff = 'ftp://ftp.ncbi.nlm.nih.gov/refseq//H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz'
req  = Request(gff)
req.add_header('Accept-Encoding', 'gzip')
response = urlopen(req)
content = gzip.decompress(response.read())
gff_file = content.splitlines()

for line in gff_file:
    line = line.decode()
    if not line.startswith('#'):
        fields = line.strip().split('\t')
        if fields[2] == "gene":
            chrom, start, end, \
               strand, extra_fields = itemgetter(0,3,4,6,-1)(fields)
            info_dict = parse_extra_fields(extra_fields)
#            line = '{chrom}\t{start}\t{end}\t{gene_name}\t'\
#                    '0\t{strand}\t{gene_type}\t{gene_id}' \
#                    .format(chrom = chrom,
#                            start = start, 
#                            end = end,
#                            strand = strand,
#                            gene_name = info_dict['gene_name'],
#                            gene_id = info_dict['gene_id'],
#                            gene_type = info_dict['gene_type'])
#            print(line, file = sys.stdout)
