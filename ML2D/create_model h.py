#only focus on z errors/x measurements
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
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(128, activation=tf.nn.relu)(input)
    output = tf.keras.layers.Dense(1, activation=tf.nn.sigmoid)(output)
    model = tf.keras.Model(inputs=input, outputs=output)

    model.compile(loss='BinaryCrossentropy', optimizer='adam', metrics=['accuracy'])


    for i in range(10000):
        print(i)
        subprocess.run(''.join(['export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:~/libtensorflow2/lib; ./simulate -s TORUS --pmin ',str(p),' -n 512 --Lmin 7 --Lmax 7 -v 0 --generate --new --fname \"/scratch/users/ladmon/ML/train_data/model_h,L=5(7),layer=3x128,epochs=10000,p=',str(p),'.csv\" --binary']), shell=True)
        df=pd.read_csv("".join(['/scratch/users/ladmon/ML/train_data/model_h,L=5(7),layer=3x128,epochs=10000,p=',str(p),'.csv']))
        X = df.values[:,0:25*2]
        Y = df.values[:,25*2:25*2+1]
        model.fit(X, Y, epochs=1,batch_size=512)

    model.save("".join(['/scratch/users/ladmon/ML/models/model_h,L=5(7),layer=3x128,epochs=10000,p=',str(p)]), save_format='tf')
    print("model_saved")

if __name__ == "__main__":
    main(sys.argv[1:])