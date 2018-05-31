#!/usr/bin/env python

import urllib.request
import urllib
from sequencing_tools.io_tools import xopen

baseurl = 'http://rnacentral.org/api/v1/rna/{rna_id}/?format=fasta'
with xopen('y_rna_id.list.gz') as yRNA:
    for id in yRNA:
        id = id.split('_')[0]
        url = baseurl.format(rna_id = id)
        html = urllib.request.urlopen(url)
        print(html.read().decode().rstrip('\n'))

