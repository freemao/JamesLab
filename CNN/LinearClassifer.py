# 7/6/18
# Chenyong Miao
# design the model
"""
the linear neural network model to segment sorghum components
"""
import numpy as np
from numpy.random import uniform
from keras.utils import to_categorical
from keras.layers import Dense, Dropout, Flatten
from keras.models import Sequential
from keras.optimizers import Adam, SGD
from keras.callbacks import EarlyStopping
from keras.callbacks import History
import pickle

def load_data(predict, target):
    '''
    load numpy array object files .npy
    '''
    try:
        predictors = np.load(predict) if predict.endswith('.npy') else np.loadtxt(predict)# 8000*243
    except IOError:
        print('numpy array object file %s does not exist.'%x)
    try:
        target = np.load(target) if target.endswith('.npy') else np.loadtxt(target) # 8000*1
    except IOError:
        print('numpy array object file %s does not exist.'%y)
    target = to_categorical(target)
    print(target.shape)
    return(predictors,target)

def get_model(predictors, num_hiden_layers, num_units, ActivFunc='relu'):
    '''
    design the model structure.
    can also try different activation functions for hidden layers: 
        ReLU(Leaky), ELU(Exponential), or Maxout. 
    For relu, caution on the dying neuron problem
    '''
    n_cols = predictors.shape[1]
    model = Sequential()
    model.add(Dense(num_units, activation=ActivFunc, input_shape=(n_cols,))) # dense layer
    for i in range(int(num_hiden_layers)-1):
        model.add(Dense(num_units, activation=ActivFunc))
    model.add(Dense(4, activation='softmax')) # output layer
    model.summary()
    return(model)

def train(x, y, lyr, unit, lr, prefix):
    predictors,target = load_data(x,y)
    lyr = int(lyr)
    unit = int(unit)
    lr = float(lr)
    model = get_model(predictors, lyr, unit)
    my_optimizer = SGD(lr=lr)
    model.compile(optimizer=my_optimizer, loss='categorical_crossentropy', metrics=['accuracy']) 
    early_stopping_monitor = EarlyStopping(patience=20)
    model_history = model.fit(predictors, target, validation_split=0.2, epochs=1000, callbacks=[early_stopping_monitor])
    model.save('%s_model_%s_%s_%s.h5'%(prefix, lyr,unit,lr)) # save model
    pickle.dump(model_history.history, open( "%s_history_%s_%s_%s.p"%(prefix, lyr,unit,lr), "wb" ) ) # save training history


import sys
if len(sys.argv)==7:
    train(*sys.argv[1:])
else:
    print('np_x, np_y, lyr, unit, lr, prefix')
