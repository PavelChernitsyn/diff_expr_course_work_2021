import numpy as np

def myFunc(x, y, v, mu, l):
    # x - время, y - длина свисающей части
    a = (1 + mu) * 9.8 * y / l - mu * 9.8
    v = v + (-9.8 * (x + 0.001) * mu + (9.8 * (x + 0.001) * y * (1 + mu))/l) - (-9.8 * x * mu + (9.8 * x * y * (1 + mu))/l)
    return a, v