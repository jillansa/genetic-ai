
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

import numpy as np

def parse_alpha_to_seq(sequence):
    output = np.arange(len(sequence))
    for i in range(0, len(sequence)):
        snippet = sequence[i]
        if snippet == 'A':
            output[i] = 0
        elif snippet == 'C':
            output[i] = 1
        elif snippet == 'T':
            output[i] = 2
        elif snippet == 'G':
            output[i] = 3
        elif snippet == 'N':
            output[i] = -1
        else:
            raise AssertionError("Cannot handle snippet: " + snippet)
    return output

def to_categorical(y, nb_classes=None):
    '''Convert class vector (integers from 0 to nb_classes)
    to binary class matrix, for use with categorical_crossentropy
    '''
    y = np.asarray(y, dtype='int32')
    if not nb_classes:
        nb_classes = np.max(y) + 1
    Y = np.zeros((len(y), nb_classes))
    for i in range(len(y)):
        if y[i] != -1:
            Y[i, y[i]] = 1.
    return Y


def do_one_hot_encoding(sequence, seq_length, f=parse_alpha_to_seq):
    X = np.zeros((sequence.shape[0], seq_length, 4))
    for idx in range(0, len(sequence)):
        X[idx] = to_categorical(f(sequence[idx]), 4)
    return X

 
def parse_protein_to_array(output, value, class_length):
    for i in range(0, class_length):        
        if value == i:
            output[i] = '1'
    return output

def do_one_hot_encoding_protein_clasif(values, class_length, f=parse_protein_to_array):
    X = np.zeros((values.shape[0], class_length))
    for idx in range(0, len(values)): 
        X[idx] = f(X[idx], values[idx], class_length)
    return X

# function to convert sequence strings into k-mer words, default size = 6 (hexamer words)
def getKmers(sequence, size):
    return [sequence[x:x+size].lower() for x in range(len(sequence) - size + 1)]