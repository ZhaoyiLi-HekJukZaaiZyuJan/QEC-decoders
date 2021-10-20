from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from keras import metrics
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils
import keras
import tensorflow as tf
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.pipeline import Pipeline
import numpy as np
import pandas as pd
import matplotlib.pyplot
import sys
import os
import subprocess

#subprocess.run('export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:~/libtensorflow2/lib; build/simulate -s TORUS --pmin 0.2  --Np 1000 -n 512 --Lmin 7 --generate --fname "train_data/data1.csv"', shell=True)

df=pd.read_csv("train_data/data1.csv")
X = df.values[:,0:49*2]
Y = df.values[:,49*2:49*2+2]

loaded_2 = tf.keras.models.load_model('models/model1')
input = tf.keras.Input(shape=(2*49,))
output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
output = tf.keras.layers.Dense(2, activation=tf.nn.sigmoid)(output)
model = tf.keras.Model(inputs=input, outputs=output)
model.compile(loss='BinaryCrossentropy', optimizer='adam', metrics=['mse',tf.keras.metrics.BinaryAccuracy()])

model.fit(X, Y, epochs=5000, batch_size=512)


# Export the model to a SavedModel
model.save('models/model1', save_format='tf')

