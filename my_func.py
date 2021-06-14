import numpy as np

def myFunc(x, y):
    dy = np.zeros((len(y)))
    dy[0] = 3*(1+x) - y[0]
    return dy