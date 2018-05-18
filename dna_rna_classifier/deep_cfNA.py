#!/usr/bin/env python

import tensorflow as tf
from matplotlib import use as mpl_use
mpl_use('Agg')
from keras.models import Sequential, model_from_json
from keras.layers import Dense, Conv1D,\
                         Flatten, MaxPool1D,  \
                         Dropout, LSTM, \
                         Bidirectional
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
from sequencing_tools.fastq_tools import reverse_complement, \
                    kmer_bag, \
                    onehot_sequence_encoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score

acceptable_chrom = list(range(1,23))
acceptable_chrom.extend(['X','Y'])
acceptable_chrom = ['chr' + str(chrom) for chrom in acceptable_chrom]
frag_size = 600
acceptable_nuc = list('ACTGN')
dna_encoder = onehot_sequence_encoder(''.join(acceptable_nuc))

def deep_model():
    model = Sequential()
    model.add(Conv1D(filters=160, 
                  kernel_size = 26,
                  strides = 1,
                  padding = 'valid',
                  input_shape = (frag_size,5),
                  activation = 'relu'))
    model.add(MaxPool1D(pool_size=50, strides=13))
    model.add(Dropout(0.2))
    model.add(Bidirectional(LSTM(64)))
    model.add(Dropout(0.5))
    model.add(Dense(25, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', 
                                optimizer='rmsprop', 
                                metrics=['binary_accuracy'])
    return model

def get_padded_seq(bed_file, fasta):
    '''
    For each record in bed file, extract the sequence, and center it
    fill up both sides to length of (seq_length) with Ns.


    Bed files need these columns:
    1. chrom
    2. start
    3. end
    4. 
    5. 
    6. strand
    7. label: (DNA or RNA)
    '''
    genome_fa = pysam.Fastafile(fasta)
    for i, line in enumerate(open(bed_file)):
        fields = line.rstrip('\n').split('\t')
        chrom, start, end, strand, label = itemgetter(0,1,2,5,-1)(fields)
        if chrom != 'chrM':
            start, end = int(start), int(end)
            seq_length = end - start
            if seq_length < frag_size:
                padding_base = frag_size - seq_length
                half_padding = int(padding_base/2)
                seq = genome_fa.fetch(chrom, start, end)
                seq = seq.upper()
                seq = half_padding * 'N' + seq + (half_padding + 1) * 'N'

            else:
                center = (end + start) / 2
                seq = genome_fa.fetch(chrom, 
                                      int(center) - int(frag_size/2), 
                                      int(center) + int(frag_size/2))
                seq = seq.upper() 
                seq = reverse_complement(seq) if strand == "-" else seq
            
            yield seq[:frag_size], label


def generate_padded_data(bed_file, fasta):
    '''
    Wrapper for generating one-hot-encoded sequences
    '''
    for i, (seq, label) in enumerate(get_padded_seq(bed_file, fasta)):
        if set(seq).issubset(acceptable_nuc):
            label = 0 if label == "DNA" else 1
            yield dna_encoder.transform(seq), label


class data_generator():
    
    def __init__(self, bed_file, fasta, batch_size=1000):
        '''
        Wrapper for generating one-hot-encoded sequences

        return batches
        '''
        self.bed = bed_file
        self.fasta = fasta
        self.batch_size = batch_size
        self.generator = get_padded_seq(self.bed, self.fasta)

    def data_gen(self):
        X, Y = [], []

        for i in range(self.batch_size):
            try:
                seq, label = next(self.generator)
            except StopIteration:
                self.generator = get_padded_seq(self.bed, self.fasta)
                seq, label = next(self.generator)

            if set(seq).issubset(acceptable_nuc):
                X.append(dna_encoder.transform(seq))
                label = 0 if label == "DNA" else 1
                Y.append(label)
        return X, Y


    def __next__(self):
        X, Y = self.data_gen()
        return np.array(X), np.array(Y)


def fetch_validation(test_bed, fa_file):
    data = [data for data in generate_padded_data(test_bed, fa_file)]
    features, labels = zip(*data)

    features = np.array(features)
    labels = np.array(labels)
    return features, labels


def load_model(prefix):
    '''
    keras load model
    '''
    json = open(prefix + '_architecture.json', 'r').read()
    model = model_from_json(json)
    model.load_weights(prefix + '_weights.h5')
    return model


def save_model(model, prefix = 'model'):
    '''
    keras save model, weight
    '''
    weight_h5 = prefix + '_weights.h5'
    model.save_weights(weight_h5)
    print('Saved weights to %s' %weight_h5)

    # Save the model architecture
    model_json = prefix + '_architecture.json'
    with open(model_json, 'w') as f:
        f.write(model.to_json())

    print('Saved model to %s' %model_json)



def training_sample(train_bed, fa_file):
    '''
    Set up keras model
    retrieve data and train
    '''
    
    model = deep_model()
    history = model.fit_generator(data_generator(train_bed, 
                                                 fa_file, 
                                                batch_size = 10000),
                                  epochs = 1,
                                  steps_per_epoch = 10)


    print('Fitted model')
    return history, model



def validation_sample(test_bed, fa_file, model):
    # validation of the model
    X_test, y_test = fetch_validation(test_bed, fa_file)
    print('Fetched test samples')
    y_pred = model.predict_classes(X_test)
    y_pred = y_pred.flatten()
    print("Precision: %1.3f" % precision_score(y_test, y_pred))
    print("Recall: %1.3f" % recall_score(y_test, y_pred))
    print("F1: %1.3f" % f1_score(y_test, y_pred))
    print("AUROC: %1.3f" % roc_auc_score(y_test, y_pred))


def plot_train(history):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key, vals in history.history.items():
        ax.plot(np.arange(len(vals)), vals, '-o',  label = key)
    ax.legend()
    fig.savefig('deep_train.png', bbox_inches='tight', transparent = True)
    return 0


def main():
    work_dir = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier'
    test_bed = work_dir + '/test.bed'
    train_bed = work_dir + '/train.bed'
    fa_file = '/stor/work/Lambowitz/ref/hg19/genome/hg19_genome.fa'
    model_prefix = work_dir + '/deef_cfNA'

    train = False

    if train:
        history, model = training_sample(train_bed, fa_file)
        save_model(model, prefix = model_prefix)

    else:
        model = load_model(model_prefix)

    validation_sample(test_bed, fa_file, model)



if __name__ == '__main__':
    main()
