#!/usr/bin/env python

import torch
from multiprocessing import cpu_count
torch.set_num_threads(cpu_count())
import torch.optim as optim
import torch.nn.functional as F
import torch.multiprocessing as mp
import numpy as np
from torch.autograd import Variable
import pyximport
pyximport.install()
from utils.bed_utils import progress, data_generator, prediction_generator
from utils.model import Deep_cfNA, calculate_metrics
import time
N_PADDED=True


def validate(test_bed, fa, model_file):
    '''
    run validation by a feature generator
    '''
    model = Deep_cfNA()
    model.load_state_dict(torch.load(model_file))
    model.eval()
    print('Loaded model: ', model_file)
    pred_Y = []
    Y = []

    for X, lines in prediction_generator(test_bed, fa, batch_size = 500, N_padded=N_PADDED):
        pred = model(X)
        pred_Y.extend(pred.reshape(-1).numpy())
        y = [line.split('\t')[-1] for line in lines]
        Y.extend(y)
    
    loss = F.binary_cross_entropy(torch.Tensor(pred_Y), torch.Tensor(y))
    calculate_metrics(Y, pred_Y, loss.item())

def deep_train(data_iterator, model, epoch=0, steps=500):
    '''
    training for one epoch
    '''
    optimizer = optim.Adam(model.parameters(), lr = 0.001)
    print('Start training epoch %i.....' %epoch)
    for step in range(steps):
        start = time.time()
        optimizer.zero_grad()

        # get data
        X, y = next(data_iterator)
        assert X.shape==(data_iterator.batch_size,5, 400) or \
                X.shape==(data_iterator.batch_size + 1,5,400), \
                X.shape

        #get prediction (forward)
        pred_y = model(X)
        pred_y = pred_y.view(-1)
        assert sum(pred_y != pred_y).item() == 0, pred_y
        loss = F.binary_cross_entropy(pred_y, y)
        
        # update gradient
        loss.backward()
        optimizer.step()
        end = time.time()

        status = '%i/%i step\tloss:%.3f\tused: %.3fs' \
                %(step, steps, loss.item(), end-start)
        progress(steps, step, epoch + 1, status)
    calculate_metrics(y, pred_y, epoch + 1, loss.item())


def train(RNA_bed, DNA_bed, fa):
    '''
    control training
    '''
    model = Deep_cfNA()
    model.train()
    model.initialize_weight()

    batch = 500
    epochs = 5
    for epoch in range(epochs):
        data_iterator = data_generator(RNA_bed, DNA_bed, fa, 
                                       batch_size=batch, 
                                       N_padded = N_PADDED,
                                       seed = epoch)
        deep_train(data_iterator, model, epoch + 1, steps = 10000)

    return model
            

def main():
    work_dir = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier'
    DNA_bed = work_dir + '/train_DNA.bed'
    RNA_bed = work_dir + '/train_RNA.bed'
    test_bed = work_dir + '/test.bed'
    fa = '/stor/work/Lambowitz/ref/hg19/genome/hg19_genome.fa'
    model_file = work_dir + '/pytorch_cfNA.pt'

    model = train(RNA_bed, DNA_bed, fa)

    if model_file:
        torch.save(model.state_dict(), model_file)
        print('Saved: ', model_file)



if __name__ == '__main__':
    main()




