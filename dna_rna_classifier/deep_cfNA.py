#!/usr/bin/env python

import tensorflow as tf
from keras.models import Sequential, model_from_json
from keras.layers import Dense, Conv1D,\
                         Flatten, MaxPool1D,  \
                         Dropout, LSTM, \
                         Bidirectional
from keras import backend as K
from keras.utils import plot_model
import pysam
import numpy as np
from operator import itemgetter
from sequencing_tools.fastq_tools import reverse_complement, \
                    kmer_bag, \
                    onehot_sequence_encoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score
from collections import defaultdict
import random

'''
Only take in fragments from regular chromosomes
'''
acceptable_chrom = list(range(1,23))
acceptable_chrom.extend(['X','Y'])
acceptable_chrom = ['chr' + str(chrom) for chrom in acceptable_chrom]
frag_size = 400
acceptable_nuc = list('ACTGN')
dna_encoder = onehot_sequence_encoder(''.join(acceptable_nuc))

def f1(y_true, y_pred):
    '''
    https://stackoverflow.com/questions/43547402/how-to-calculate-f1-macro-in-keras?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    '''
    def recall(y_true, y_pred):
        """Recall metric.

        Only computes a batch-wise average of recall.

        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

    def precision(y_true, y_pred):
        """Precision metric.

        Only computes a batch-wise average of precision.

        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

def deep_model():
    '''
    DanQ model 
    https://github.com/uci-cbcl/DanQ/blob/master/DanQ-JASPAR_train.py
    '''
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
                metrics=[f1,'binary_accuracy'])
    plot_model(model, to_file='model.png', show_shapes=True)
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
            label = 1 if label == "DNA" else 0
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
        '''
        Populate reponse vector and feature array with desired batch size
        '''
        X, Y = [], []

        label_counter = defaultdict(int)
        for i in range(self.batch_size):
            try:
                seq, label = next(self.generator)
            except StopIteration:
                self.generator = get_padded_seq(self.bed, self.fasta)
                seq, label = next(self.generator)

            if set(seq).issubset(acceptable_nuc):
                if label_counter[label] <= self.batch_size/2 and random.random() >= 0.5:
                    X.append(dna_encoder.transform(seq))
                    label = 1 if label == "DNA" else 0
                    Y.append(label)
                    label_counter[label] += 1
        return X, Y


    def __next__(self):
        '''
        generator for Keras fit_generator
        '''
        X, Y = self.data_gen()
        return np.array(X), np.array(Y)


def fetch_validation(test_bed, fa_file, batch_size = 0):
    '''
    fetch sequences from test bed file and return feature arrays and test label
    '''
    label_count = defaultdict(int)
    if batch_size > 0:
        features = []
        labels = []
        for seq, label in generate_padded_data(test_bed, fa_file):
            if label_count[label] < batch_size/2:
                features.append(seq)
                labels.append(label)
                label_count[label] += 1
    
    else:
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
                                                batch_size = 500),
                                  epochs = 20,
                                  steps_per_epoch = 1000)


    print('Fitted model')
    return history, model



def validation_sample(test_bed, fa_file, model):
    '''
    Test model on unseen data
    '''
    # validation of the model
    X_test, y_test = fetch_validation(test_bed, fa_file)
    print('Fetched test samples')

    # make prediction
    y_pred_prob = model.predict_proba(X_test)
    y_pred_prob = y_pred_prob.flatten()

    y_pred_class = model.predict_classes(X_test)
    y_pred_class = y_pred_class.flatten()

    # evaluation
    print("Precision: %1.3f" % precision_score(y_test, y_pred_class))
    print("Recall: %1.3f" % recall_score(y_test, y_pred_class))
    print("F1: %1.3f" % f1_score(y_test, y_pred_class))
    print("AUROC: %1.3f" % roc_auc_score(y_test, y_pred_prob))


def plot_train(history):
    from matplotlib import use as mpl_use
    mpl_use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt

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

    train = True

    if train:
        history, model = training_sample(train_bed, fa_file)
        save_model(model, prefix = model_prefix)
        plot_train(history)

    else:
        model = load_model(model_prefix)

    validation_sample(test_bed, fa_file, model)



if __name__ == '__main__':
    main()
