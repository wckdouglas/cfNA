#!/usr/bin/env python

import pymysql
import pandas as pd
from urllib.request import Request, urlopen
import gzip
from operator import itemgetter


def parse_extra_fields(extra_field):
    info_fields = extra_field.split(';')
    info_dict = {}
    for info in info_fields:
        row_fields = info.strip().split('=')
        info_dict[row_fields[0]] = row_fields[1].strip('')
    return info_dict


def fetch_gtf():
    gff = 'ftp://ftp.ncbi.nlm.nih.gov/refseq//H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz'
    req  = Request(gff)
    req.add_header('Accept-Encoding', 'gzip')
    response = urlopen(req)
    content = gzip.decompress(response.read())
    gff_file = content.splitlines()

    rows = []
    for line in gff_file:
        line = line.decode()
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if fields[2] == "gene":
                start, end, \
                    strand, extra_fields = itemgetter(3,4,6,-1)(fields)
                info_dict = parse_extra_fields(extra_fields)

                start = int(start) - 1
                end = int(end) 
                rows.append((start, end, strand, info_dict['gene_biotype'], info_dict['gene']))    
    
    return pd.DataFrame(rows, columns = ['start','end', 'strand','gene_type','gene_name'])
    

def fetch_bed():
    conn = pymysql.connect(host = 'genome-mysql.soe.ucsc.edu',
                    user ='genome',
                    db = 'hg19',
                    port = 3306)

    return pd.read_sql_query('SELECT rs.chrom, rs.txstart, rs.txEnd, rs.name, rs.score, rs.strand, rs.name2 '\
                    'FROM ncbiRefSeq as rs', 
                    conn)  \
            .rename(columns = {'name2':'gene_name',
                                'txstart':'start',
                                'txEnd':'end' })


def main(refseq_bed):
    gtf = fetch_gtf()
    bed = fetch_bed()
    bed.merge(gtf, how = 'inner', 
            on = ['start','end','strand','gene_name']) \
        .pipe(lambda d: d[['chrom','start','end','gene_name',
                        'score','strand','gene_type','name']]) \
        .to_csv(refseq_bed, 
                sep ='\t',
                index=False,
                header=False)


if __name__ == '__main__':
    main('/stor/work/Lambowitz/ref/hg19_ref/genes/hg19_refseq.bed')

