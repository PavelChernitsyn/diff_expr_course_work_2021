import numpy as np


def myFunc(x, y, v, mu, l, h):
    # x - время, y - длина свисающей части
    a = (1 + mu) * 9.8 * y / l - mu * 9.8
    v = v + (-9.8 * (x + h) * mu + (9.8 * (x + h) * y * (1 + mu))/l) - (-9.8 * x * mu + (9.8 * x * y * (1 + mu))/l)
    return a, v