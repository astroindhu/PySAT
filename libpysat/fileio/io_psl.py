import os
import numpy as np

def fold2dict(path):
    '''Covert the folder of text files to python dictionary'''
    data = {}
    for f in os.listdir(path):
        if f.endswith("txt"):
            data[f] = np.loadtxt(path + f)
    return data