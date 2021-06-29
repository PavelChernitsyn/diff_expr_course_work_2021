import numpy as np

def myFunc(y, env_temp):
    dy = np.zeros((len(y)))
    # dy[0] = 3*(1+x) - y[0]
    dy[0] = 1/25 * np.log(2) * (env_temp - y)
    return dy