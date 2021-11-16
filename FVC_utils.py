import pickle
from os import path, makedirs
import time
import pandas as pd
import numpy as np

def printx(output):
    print(output, flush = True)


def save_to_pickle(var, filename='temp', root_path='./'):
    '''
    Save variable to pickle file defaultly in current dir
    '''
    #fn = path.join(root_path, filename)
    fn = filename
    if path.exists(path.split(fn)[0]) == False:
        makedirs(path.split(fn)[0])
    #printx("saving {} into {}".format(filename, root_path))
    printx("saving into {}".format(filename))
    f = open(fn, 'wb')
    pickle.dump(var, f)
    f.close()
    return var


def load_from_pickle(filename='temp', root_path='./'):
    '''
    Load pickle file defaultly from current dir
    '''
    #printx("loading {} from {}".format(filename, root_path))
    #f = open(path.join(root_path, filename), 'rb')
    printx("loading {}".format(filename))
    f = open(filename, 'rb')
    var = pickle.load(f)
    f.close()
    return var


def readData(file_path1, file_path2):
    # Read file and combine
    data1 = pd.read_csv(file_path1, header=None, sep=' ', low_memory=False)
    data2 = pd.read_csv(file_path2, header=None, sep=' ', low_memory=False)
    
    combined = pd.concat([data1, data2])
    # Return X, y
    return combined.iloc[:, 6:].values, combined.iloc[:, 0].values
