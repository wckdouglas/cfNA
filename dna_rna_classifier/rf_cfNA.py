#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import pysam
from operator import itemgetter
from sequencing_tools.fastq_tools import reverse_complement, kmer_bag, onehot_sequence_encoder
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score
from sklearn.preprocessing import Binarizer, label_binarize
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np


test_bed = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier/no_snc_test.bed'
fa_file = '/stor/work/Lambowitz/ref/hg19/genome/hg19_genome.fa'


acceptable_chrom = list(range(1,23))
acceptable_chrom.extend(['X','Y'])
acceptable_chrom = ['chr' + str(chrom) for chrom in acceptable_chrom]

def get_seq(bed_file, fasta):
    genome_fa = pysam.Fastafile(fasta)
    for line in open(bed_file):
        fields = line.rstrip('\n').split('\t')
        chrom, start, end, strand, label = itemgetter(0,1,2,5, -1)(fields)
        if chrom in acceptable_chrom:
            seq = genome_fa.fetch(chrom, int(start), int(end))
            seq = reverse_complement(seq) if seq == "-" else seq
            seq = seq.upper() 
            yield seq, label


def extract_data(bed_file, fasta):
    label_counter = defaultdict(int)
    for i, (seq, label) in enumerate(get_seq(bed_file, fasta)):
        if 'N' not in seq:
            kmers = kmer_bag(seq, k_start =  2, k_end = 6)
            outid = 'read%i' %(i+1)
            label_counter[label] += 1
            if label_counter[label] < 40000:
                for kmer, kmer_count in kmers.items():
                    yield(outid, kmer, kmer_count, label)


data_generator = extract_data(test_bed, fa_file)
kdf = pd.DataFrame(data_generator,
                  columns = ['read_id','kmer','kmer_count','label'])
print('Fetched data')


sdf = kdf\
    .pipe(pd.pivot_table, 
        index = ['read_id','label'],
        columns = 'kmer',
        values = 'kmer_count',
        fill_value = 0) \
    .reset_index()

X = sdf.filter(regex = '[ACTG]')
Y = sdf.label
print('%i samples with labels: ' %(Y.shape[0]), Y.unique())

X_train, X_test, y_train, y_test = train_test_split(X, Y, 
                                                    test_size=0.2, 
                                                    random_state=0)
rf = RandomForestClassifier(oob_score=False)
rf_cv = GridSearchCV(rf, 
                     n_jobs = 24,
                     param_grid = {'n_estimators': np.arange(10,30)})
rf_cv.fit(X_train, y_train)
print('Trained model')


y_pred = rf_cv.best_estimator_.predict_proba(X_test)[:, 1].round()
y_test_code = y_test.astype('category').cat.codes

print("Precision: %1.3f" % precision_score(y_test_code, y_pred))
print("Recall: %1.3f" % recall_score(y_test_code, y_pred))
print("F1: %1.3f" % f1_score(y_test_code, y_pred))
print("AUROC: %1.3f" % roc_auc_score(y_test_code, y_pred))


fig = plt.figure()
ax = fig.add_subplot(111)
sns.barplot(data = pd.DataFrame({'features':X.columns,
                      'importance': rf_cv.best_estimator_.feature_importances_
                     }) \
                .sort_values('importance', ascending=False)\
                .head(30),
            x = 'features',
            y = 'importance', 
            color = 'steelblue',
            ax = ax)
xt = ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
sns.despine()
fig.savefig('rf_factor.png')
