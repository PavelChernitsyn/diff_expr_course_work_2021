import numpy as np


def myFunc(x, mu, l):
    # x - длина свисающей части
    a = (1 + mu) * 9.8 * x[0] / l - mu * 9.8
    # v = v + (-9.8 * (x + h) * mu + (9.8 * (x + h) * y * (1 + mu))/l) - (-9.8 * x * mu + (9.8 * x * y * (1 + mu))/l)
    # v += a*h

    return x[-1], a 