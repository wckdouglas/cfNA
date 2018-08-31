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
                                 out_channels=160, #(10 neurons/filters)
                                kernel_size=26, #(scanning 26 nucleotides at a time)
                                stride=1) #(moving 1 nucleotide at a time) 
        self.max_pool = nn.MaxPool1d(kernel_size = 50, stride = 13)
        """
        self.LSTM = nn.LSTM(input_size=160,
                          hidden_size=64, 
                          batch_first = True,
                          bidirectional=True) #(64 neurons)
        self.linear1 = nn.Linear(128, 50) #(128 input from 2x64 LSTM)
        """
        self.linear0 = nn.Linear(4160, 160)
        self.linear1 = nn.Linear(160, 25)
        self.linear2 = nn.Linear(25, 10) 
        self.linear3 = nn.Linear(10, 1)


    def initialize_weight(self):
        for name, param in self.named_parameters():
            if 'bias' in name:
                nn.init.constant_(param, 0.0)
            elif 'weight' in name:
                nn.init.xavier_normal_(param)


    def forward(self, x):
        '''
        get prediction, needed for pytorch.nn.Module
        input: (batch_size, num_feature, seq_length)
        '''

        y = self.conv_1d(x)
        y = F.relu(y)
        y = self.max_pool(y)
        y = F.dropout(y, p = 0.2)

        """
        '''
        #change input to (batch_size, seq_length, num_feature) for LSTM layer
        #'''
        y = y.transpose(1,2)
        y, (hidden, cell) = self.LSTM(y) # bidirectional LSTM concat
        y = y[:,-1]  # last time stamp from LSTM layer
        y = F.dropout(y, p = 0.5)
        """
        batch_size = y.size(0)
        y = y.view(batch_size, -1)
        y = self.linear0(y)
        y = F.dropout(y, p = 0.2)
        y = self.linear1(y)
        #y = self.linear2(y)
        y = F.dropout(y, p = 0.2)
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
          ' Loss: ', round(loss,3), 
          'Accuracy: ', round(accuracy, 3), 
          ' F1: ', round(f1,3), 
          ' Precision: ', round(precision,3), 
          ' Recall: ', round(recall,3))
