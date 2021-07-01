import numpy as np


def myFunc(x, v, mu, l):
    # x - длина свисающей части
    a = (1 + mu) * 9.8 * x / l - mu * 9.8
    # v = v + (-9.8 * (x + h) * mu + (9.8 * (x + h) * y * (1 + mu))/l) - (-9.8 * x * mu + (9.8 * x * y * (1 + mu))/l)
    # v += a*h
    return a, v