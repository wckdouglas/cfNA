#!/usr/bin/env python

import torch
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
from torch.autograd import Variable
import pyximport
pyximport.install()
from utils.bed_utils import progress, data_generator, prediction_generator
from utils.model import Deep_cfNA, calculate_metrics
import torch.multiprocessing as mp
torch.set_num_threads(24)


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

    for X, lines in prediction_generator(test_bed, fa, batch_size = 500):
        pred = model(X)
        pred_Y.extend(pred.reshape(-1).numpy())
        y = [line.split('\t')[-1] for line in lines]
        Y.extend(y)
    
    loss = F.binary_cross_entropy(torch.Tensor(pred_Y), torch.Tensor(y))
    calculate_metrics(Y, pred_Y, loss.item())


def deep_train(rank, data_iterator, model, epochs = 5, steps=500):
    optimizer = optim.RMSprop(model.parameters(), lr = 0.001)
    print('[CPU %i] Start training.....' %(rank + 1))

    for epoch in range(epochs):
        for step in range(steps):
            X, y = next(data_iterator)
            assert X.shape==(data_iterator.batch_size,5, 400) or \
                    X.shape==(data_iterator.batch_size + 1,5,400), \
                    X.shape
            pred_y = model(Variable(X))
            pred_y = pred_y.reshape(-1)
            loss = F.binary_cross_entropy(pred_y, y)
            
            # update gradient
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            progress(steps, step, epoch + 1, epochs)
        calculate_metrics(y, pred_y, epoch + 1, loss.item())



def train(rank, RNA_bed, DNA_bed, fa, model):
    batch = 200
    data_iterator = data_generator(RNA_bed, DNA_bed, fa, 
                                   batch_size=batch, N_padded = False)
    deep_train(rank, data_iterator, model, epochs = 5)

def train_wrapper(RNA_bed, DNA_bed, fa, model_file):
    model = Deep_cfNA()
    model.share_memory()

    processes = []
    for rank in range(2):
        p = mp.Process(target=train, args=(rank, RNA_bed, DNA_bed, fa, model))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()


    if model_file:
        torch.save(model.state_dict(), model_file)
        print('Saved: ', model_file)
            

def main():
    work_dir = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier'
    DNA_bed = work_dir + '/train_DNA.bed'
    RNA_bed = work_dir + '/train_RNA.bed'
    test_bed = work_dir + '/test.bed'
    fa = '/stor/work/Lambowitz/ref/hg19/genome/hg19_genome.fa'
    model_file = work_dir + '/pytorch_cfNA.pt'

    train_wrapper(RNA_bed, DNA_bed, fa, model_file)



if __name__ == '__main__':
    main()




