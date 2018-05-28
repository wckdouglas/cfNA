import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score

class Deep_cfNA(nn.Module):
    '''
    A model with DanQ architecture, but in a smaller scale

    Keras implementation: https://github.com/wckdouglas/deep_cfNA/blob/master/deep_cfNA/deep_cfNA_model.py
    '''
    def __init__(self):
        super(Deep_cfNA, self).__init__()
        self.conv_1d = nn.Conv1d(in_channels=5, #(ACTGN)
                                 out_channels=160, #(160 neurons/filters)
                                kernel_size=26, #(scanning 26 nucleotides at a time)
                                stride=1) #(moving 1 nucleotide at a time) 
        self.LSTM = nn.LSTM(input_size=26,hidden_size=64, bidirectional=True) #(64 neurons)
        self.linear1 = nn.Linear(128, 50) #(128 input from 2x64 LSTM)
        self.linear2 = nn.Linear(50, 25) 
        self.linear3 = nn.Linear(25, 1)

        self.linear1A = nn.Linear(26, 128)
        self.linear1B = nn.Linear(128, 1)

    
    def forward(self, x):
        y = self.conv_1d(x)
        y = F.relu(y)
        y = F.max_pool1d(y, kernel_size=50, stride=13)
        y = F.dropout(y, p = 0.2)
        y, (hidden, cell) = self.LSTM(y)
        y = y[:,-1]
        y = F.dropout(y, p = 0.5)
        y = self.linear1(y)
        y = self.linear2(y)
        y = self.linear3(y)
        y = F.sigmoid(y)
        return y

def calculate_metrics(y, pred_y, epoch, loss):
    '''
    Output some metrics
    '''
    pred_label = pred_y.detach().numpy().round()
    accuracy = accuracy_score(y, pred_label)
    precision = precision_score(y, pred_label)
    recall = recall_score(y, pred_label)
    f1 = f1_score(y, pred_label)
    print('\nMini-batch: ', epoch, 
          ' Loss: ', loss, 
          'Accuracy: ', accuracy, 
          ' F1: ', f1, 
          ' Precision: ', precision, 
          ' Recall: ', recall)
