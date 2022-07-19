from keras.utils import np_utils
import tensorflow as tf
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import argparse
import subprocess
import sys

def main(argv):
    parser = argparse.ArgumentParser(description='Simulation parser')
    parser.add_argument('-p','--probability',default = "0.001")
    parser.parse_args()
    args = parser.parse_args()
    p = args.probability
    print("probability:", p, "\n")

    input = tf.keras.Input(shape=(2*25,))
    output = tf.keras.layers.Dense(64, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(64, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(16, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(4, activation=tf.nn.softmax)(output)
    model = tf.keras.Model(inputs=input, outputs=output)
    myopt = tf.keras.optimizers.Adam()
    model.compile(loss='categorical_crossentropy', optimizer=myopt, metrics=['accuracy'])


    for i in range(10000):
        print(i)
        subprocess.run(''.join(['export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:~/libtensorflow2/lib; ./simulate -s TORUS --pmin ',str(p),' -n 512 --Lmin 7 --Lmax 7 -v 0 --generate --new --fname \"/scratch/users/ladmon/ML/train_data/model_shuttle,L=5(7),layer=3x128,epochs=10000,p=',str(p),'.csv\"']), shell=True)
        df=pd.read_csv("".join(['/scratch/users/ladmon/ML/train_data/model_shuttle,L=5(7),layer=3x128,epochs=10000,p=',str(p),'.csv']))
        X = df.values[:,0:25*2]
        Y = df.values[:,25*2:25*2+1]
        # one hot encoding
        # encode class values as integers
        encoder = LabelEncoder()
        encoder.fit(Y)
        encoded_Y = encoder.transform(Y)
        # convert integers to dummy variables (i.e. one hot encoded)
        dummy_y = np_utils.to_categorical(encoded_Y, num_classes=4)
        model.fit(X, dummy_y, epochs=1,batch_size=512)


    #=======================================================================================

    #df=pd.read_csv("train_data/four_ways_data_0.01_w5.csv")
    #X = df.values[:,0:25*2]
    #Y = df.values[:,25*2:25*2+2]
    #
    #encoder = LabelEncoder()
    #encoder.fit(Y)
    #encoded_Y = encoder.transform(Y)
    ## convert integers to dummy variables (i.e. one hot encoded)
    #dummy_y = np_utils.to_categorical(encoded_Y)
    #


    #=======================================================================================
    #df=pd.read_csv("train_data/4_ways_data_0.08.csv")
    #X = df.values[:,0:49*2]
    #Y = df.values[:,49*2:49*2+1]
    #encoder = LabelEncoder()
    #encoder.fit(Y)
    #encoded_Y = encoder.transform(Y)
    ## convert integers to dummy variables (i.e. one hot encoded)
    #dummy_y = np_utils.to_categorical(encoded_Y)
    #
    #input = tf.keras.Input(shape=(2*49,))
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(4, activation=tf.nn.softmax)(output)
    #model2 = tf.keras.Model(inputs=input, outputs=output)
    #model2.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #model2.fit(X, dummy_y, epochs=200, batch_size=512)
    #=======================================================================================
    ##loaded_2 = tf.keras.models.load_model('models/model1')
    #input = tf.keras.Input(shape=(2*49,))
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(512, activation=tf.nn.relu)(input)
    #output = tf.keras.layers.Dense(2, activation=tf.nn.softmax)(output)
    #model = tf.keras.Model(inputs=input, outputs=output)
    #model.compile(loss='BinaryCrossentropy', optimizer='adam', metrics=['mse',tf.keras.metrics.BinaryAccuracy()])
    #
    #model.fit(X, Y, epochs=1000, batch_size=512)

    #=======================================================================================
    #Export the model to a SavedModel
    model.save("".join(['/scratch/users/ladmon/ML/models/model_shuttle,L=5(7),layer=3x128,epochs=10000,p=',str(p)]), save_format='tf')
    print("model_saved")

if __name__ == "__main__":
    main(sys.argv[1:])