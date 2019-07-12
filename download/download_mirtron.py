from bs4 import BeautifulSoup
import pandas as pd
import urllib
import os

mirtron_link = 'http://compgen.cshl.edu/mirna/mam_mirtrons/hg19_candidate/'
url = urllib.request.urlopen(mirtron_link)
bs = BeautifulSoup(url)
tab = bs.find('table')
pd.read_html(tab.decode())[0] \
    .to_csv(os.environ['REF'] + '/hg19_ref/genes/mirtron.csv', index=False)
