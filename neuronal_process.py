
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

import pandas as pd
import numpy as np

import tensorflow as tf
from keras.models import Model, Sequential
from keras.layers import Input, SimpleRNN, Dense, Flatten, Embedding, LSTM, Bidirectional, Dropout
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
from keras.optimizers import SGD, Adam, Adadelta, RMSprop

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
from sklearn.metrics import confusion_matrix
from sklearn.naive_bayes import MultinomialNB
from sklearn.feature_extraction.text import CountVectorizer

from neuronal_helper import *

def createClasiffierLSTM(seq_data):

    seq_data['words'] = seq_data.apply(lambda x: getKmers(x['SEQUENCE'],6), axis=1)
    seq_data = seq_data.drop('SEQUENCE', axis=1)

    #seq_data.head()
    seq_texts = list(seq_data['words'])
    for item in range(len(seq_texts)):
        seq_texts[item] = ' '.join(seq_texts[item])
    labels = seq_data.iloc[:, 2].values # PROTEIN_ID

    tokenizer = Tokenizer() 
    tokenizer.fit_on_texts(seq_texts) # tokenizer the word in each secuence
    encoded_docs = tokenizer.texts_to_sequences(seq_texts) # Transform ¿unique? each token in a integer value
    max_length = max([len(s) for s in encoded_docs]) # 18921 max langth of all secuences
    X = pad_sequences(encoded_docs, maxlen = max_length, padding = 'post') # the context is determinate in less 100 nucleotid

    X_train,X_test,y_train,y_test = train_test_split(X,labels,
                                                    test_size=0.20,random_state=42)
    vocab_size = len(tokenizer.word_index) + 1
    print(X_train.shape)
    print(X_test.shape)

    # Keras ofrece una capa de incrustación que se puede utilizar para redes neuronales en datos de texto. 
    # Requiere que los datos de entrada estén codificados en números enteros, 
    # de modo que cada palabra esté representada por un número entero único.
    #n_timesteps, n_features, n_outputs = X_train.shape[0], X_train.shape[1], y_train[0].shape

    model = Sequential()
    model.add(Embedding(vocab_size, 32))
    model.add(LSTM(32))
    model.add(Dense(1, activation='sigmoid'))

    model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['acc'])
    #history_rnn = model.fit(texts_train, y_train, epochs=10, batch_size=60, validation_split=0.2)

    epochs = 5
    #model.compile(loss='mean_absolute_error',optimizer='adam',metrics=['accuracy'])
    #binary_crossentropy
    stringlist = []
    model.summary(print_fn=lambda x: stringlist.append(x))
    short_model_summary = "\n".join(stringlist)
    print(short_model_summary)

    history = model.fit(X_train, y_train, 
                        epochs = epochs, verbose = 1, validation_split = 0.2, 
                        batch_size = 32)

    pred = model.predict_classes(X_test)
    acc = model.evaluate(X_test, y_test)

    print("Test loss is {0:.2f} accuracy is {1:.2f} % ".format(acc[0],acc[1]*100))
    print(confusion_matrix(pred, y_test))

    statistics = "" 
    statistics += "<div>Confusion matrix</div>"  
    dfStats = pd.DataFrame(confusion_matrix(pred, y_test))
    statistics += "<div>" + dfStats.to_html() + "</div>"
    statistics += "<div> accuracy = "+str(acc[1]*100) + " % </div>"

    return tokenizer, model, statistics, short_model_summary

def clasiffierLSTM(cbp, tokenizer, model):

    seq_texts = list(getKmers(cbp.querry_seq,6))
    seq_texts = ' '.join(seq_texts)

    seq = list()
    seq.append(seq_texts)

    # encode document
    X_seq = tokenizer.texts_to_sequences(seq)

    y_pred = model.predict_classes(X_seq)
    print("Predictions x_test: ", str(y_pred[0]))

    if y_pred[0]==1:
        return "Identificate pathogen: SARS-CoV-2 | Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1(complete genome)"
    else:
       return "No identificate pathogen SARS-CoV-2"


def createClassifierMNB(seq_data):
       
    seq_data['words'] = seq_data.apply(lambda x: getKmers(x['SEQUENCE'],6), axis=1)
    seq_data = seq_data.drop('SEQUENCE', axis=1)

    seq_data.head()
    seq_texts = list(seq_data['words'])
    for item in range(len(seq_texts)):
        seq_texts[item] = ' '.join(seq_texts[item])
    y_data = seq_data.iloc[:, 2].values # PROTEIN_ID

    # Now we will apply the BAG of WORDS using CountVectorizer using NLP
    # Creating the Bag of Words model using CountVectorizer()
    # This is equivalent to k-mer counting
    # The n-gram size of 4 was previously determined by testing
    cv = CountVectorizer(ngram_range=(4,4))
    X = cv.fit_transform(seq_texts)

    print(X.shape)
    
    # Splitting the human dataset into the training set and test set
    X_train, X_test, y_train, y_test = train_test_split(X, 
                                                        y_data, 
                                                        test_size = 0.20, 
                                                        random_state=42)

    print(X_train.shape)
    print(X_test.shape)
    #(3504, 232414)
    #(876, 232414)

    ### Multinomial Naive Bayes Classifier ###
    # The alpha parameter was determined by grid search previously
    classifier = MultinomialNB(alpha=0.1)
    classifier.fit(X_train, y_train)

    MultinomialNB(alpha=0.1, class_prior=None, fit_prior=True)

    y_pred = classifier.predict(X_test)
    print("Predictions x_test: ", y_pred[0:100])

    print("Confusion matrix\n")
    print(pd.crosstab(pd.Series(y_test, name='Actual'), pd.Series(y_pred, name='Predicted')))

    def get_metrics(y_test, y_predicted):
        accuracy = accuracy_score(y_test, y_predicted)
        precision = precision_score(y_test, y_predicted, average='weighted')
        recall = recall_score(y_test, y_predicted, average='weighted')
        f1 = f1_score(y_test, y_predicted, average='weighted')
        return accuracy, precision, recall, f1

    accuracy, precision, recall, f1 = get_metrics(y_test, y_pred)
    print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))

    statistics = "" 
    statistics += "<div>Confusion matrix</div>"
    statistics += "<div>" + pd.crosstab(pd.Series(y_test, name='Actual'), 
        pd.Series(y_pred, name='Predicted')).to_html() + "</div>"
    statistics += "<div> accuracy = "+str(accuracy*100) + " % </div>" + "<div> precision = "+str(precision*100) + " % </div>" + "<div> recall = "+str(recall*100) + " % </div>" + "<div>"

    return cv, classifier, statistics

def clasiffierMNB(cbp, cv, classifier):
 
    seq_texts = list(getKmers(cbp.querry_seq,6))
    seq_texts = ' '.join(seq_texts)

    seq = list()
    seq.append(seq_texts)

    # encode document
    X_seq = cv.transform(seq)

    print(X_seq.shape)
    y_pred = classifier.predict(X_seq)
    print("Predictions x_test: ", str(y_pred[0]))

    return str(y_pred[0])